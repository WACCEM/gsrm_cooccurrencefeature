#!/usr/bin/env python
"""
Run Analyses 1 and 2 of the COF pipeline for any set of sources, in dependency order, on one node.

For every source the steps are

    s1 (make_mcs_swath_masks) -> s2 (combine_tracking_masks) -> s3 (make_cooccurrence_masks) -> monthly (Analysis 1)
    s1 -> thresholds (calc_extreme_precip_thresholds), then s3 + thresholds -> attribution (Analysis 2)

(the thresholds of IMERG read their own input store and need no earlier step). A step starts only when the steps it needs have
finished with exit status 0; ready steps start longest-remaining-chain first, as many at a time as fit in the CPU-slot budget
(--max-slots), one every --stagger-sec seconds and only while the node has --min-free-gb of memory available. A failed step skips
the later steps of its source only; the other sources continue. Sources, commands, slots and estimated durations are in
config/config_pipeline.yaml.

Everything is written under --data-root (the scripts read COF_DATA_ROOT, see src/cof_paths.py), so a test run never touches
production. Each successful step leaves a marker <data-root>/pipeline_state/<source>/<step>.json; --resume skips steps whose marker
says success (and whose upstream steps did not run again), and an output that exists without a marker of this runner is never
overwritten unless --force is given.

Examples
  # everything for all six sources into a test area (dry run first: graph, commands, estimated schedule)
  python scripts/run_cof_pipeline.py --data-root /pscratch/sd/w/wcmca1/hackathon/tmp/round2 --dry-run
  python scripts/run_cof_pipeline.py --data-root /pscratch/sd/w/wcmca1/hackathon/tmp/round2
  # one source, or a few, only Analysis 2, continue an interrupted run
  python scripts/run_cof_pipeline.py --data-root DIR --sources scream icon
  python scripts/run_cof_pipeline.py --data-root DIR --sources imerg --analysis 2
  python scripts/run_cof_pipeline.py --data-root DIR --resume
  # only some steps (the steps they need must exist already)
  python scripts/run_cof_pipeline.py --data-root DIR --sources um --steps s3 monthly

Exit status: 0 all steps done, 1 a step failed, was skipped because of a failure or the run was interrupted, 2 usage error,
refused overwrite or failed preflight.

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import argparse
import copy
import datetime
import json
import os
import shlex
import signal
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import yaml

try:
    import psutil
except ImportError:  # the sampler and the memory gate then fall back to /proc/meminfo
    psutil = None

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"
CONFIG_DIR = REPO / "config"
sys.path.insert(0, str(SCRIPTS))
from check_zarr_store import scan_store, summarize  # noqa: E402

STEP_ORDER = ["s1", "s2", "s3", "monthly", "thresholds", "attribution"]
ANALYSIS_STEPS = {"1": ["s1", "s2", "s3", "monthly"], "2": ["s1", "s2", "s3", "thresholds", "attribution"], "both": STEP_ORDER}
STEP_ALIASES = {"1": "s1", "step1": "s1", "2": "s2", "step2": "s2", "3": "s3", "step3": "s3", "thr": "thresholds",
                "extreme": "thresholds", "attr": "attribution"}
STEP_DESCRIPTIONS = {"s1": "Step 1 swath masks", "s2": "Step 2 combined masks", "s3": "Step 3 COF masks", "monthly": "monthly rain map",
                     "thresholds": "extreme thresholds", "attribution": "extreme attribution"}


# ---------------------------------------------------------------------------------------------------------------------------
# Registry, tasks
# ---------------------------------------------------------------------------------------------------------------------------
def load_registry(path):
    """Read config_pipeline.yaml; returns (defaults, {source: source dict with 'step_cfg' per step})."""
    with open(path) as f:
        reg = yaml.safe_load(f)
    defaults, steps, sources = reg["defaults"], reg["steps"], {}
    for name, src in reg["sources"].items():
        src = dict(src)
        src["name"] = name
        cfg = {}
        for step in STEP_ORDER:
            d = dict(steps[step])
            d.update((src.get("steps") or {}).get(step, {}))
            d.setdefault("timeout_min", defaults["timeout_min"][step])
            cfg[step] = d
        src["step_cfg"] = cfg
        sources[name] = src
    return defaults, sources


def resolve_sources(requested, sources):
    """Canonical source names for the requested names or aliases (case-insensitive); all sources when none are requested."""
    if not requested:
        return list(sources)
    lookup = {}
    for name, src in sources.items():
        for alias in [name] + list(src.get("aliases", [])) + [src["config_key"], src["source_name"]]:
            lookup[alias.lower()] = name
    out = []
    for r in requested:
        if r.lower() not in lookup:
            raise SystemExit(f"Unknown source '{r}'. Known: " + ", ".join(f"{n} ({', '.join(s.get('aliases', []))})" for n, s in sources.items()))
        if lookup[r.lower()] not in out:
            out.append(lookup[r.lower()])
    return out


def normalize_steps(names):
    out = []
    for n in names:
        n = STEP_ALIASES.get(n.lower(), n.lower())
        if n not in STEP_ORDER:
            raise SystemExit(f"Unknown step '{n}'. Steps: {', '.join(STEP_ORDER)}")
        out.append(n)
    return out


class Task:
    """One step of one source: what to run, what it needs, what it makes."""

    def __init__(self, source, step, argv, env, slots, est_min, timeout_min, outputs, deps):
        self.source, self.step, self.argv, self.env = source, step, argv, env
        self.slots, self.est_min, self.timeout_min = slots, est_min, timeout_min
        self.outputs, self.deps = outputs, deps          # deps: list of (source, step)
        self.key = (source, step)
        self.will_run, self.note = True, ""
        # filled in while running
        self.start = self.end = None
        self.rc = None
        self.status = "pending"                          # pending | running | ok | failed | skipped | interrupted | done
        self.peak_rss_gb = 0.0
        self.proc = None
        self.log_path = None

    @property
    def label(self):
        return f"{self.source}/{self.step}"


def step_outputs(src, step, root, pcts):
    """Output paths of one step of one source (they decide what counts as done and what would be overwritten)."""
    sname = src["source_name"]
    return {"s1": [f"{root}mcs_masks/{sname}_mcs_masks_hp8.zarr"],
            "s2": [f"{root}all_masks/{sname}_allmasks_hp8_v1.zarr"],
            "s3": [f"{root}cof_masks/{sname}_cofmasks_hp8_v1.zarr"],
            "monthly": [f"{root}cof_masks/stats/monthly/{sname}_monthly_rainmap_cof_hp8_v1.nc"],
            "thresholds": [f"{root}extreme_precip/{sname}_precip_percentiles_6h_hp8_v1.nc"],
            "attribution": [f"{root}extreme_precip/{sname}_stormtype_spatial_p{p}.nc" for p in pcts]}[step]


def _thresholds_input(src, root):
    th = src.get("thresholds", {})
    if th.get("input_zarr"):
        return th["input_zarr"], th.get("input_var")
    return f"{root}mcs_masks/{src['source_name']}_mcs_masks_hp8.zarr", "tot_pr"


def build_tasks(sources_sel, steps_sel, sources, defaults, root, python, era5_zarr=None, extra_args=None):
    """
    Tasks for the selected sources and steps with their commands, outputs and dependencies (dependencies on steps that are not
    selected are kept: they are checked against what already exists). extra_args = {step: [arguments]} are appended to the
    command of that step (a later occurrence of an option overrides the earlier one), for tests and special runs.
    """
    cfg_sources = str(CONFIG_DIR / "config_sources.yaml")
    pcts = defaults["percentiles"]
    env_base = {k: str(v) for k, v in defaults.get("env", {}).items()}
    tasks = {}
    for name in sources_sel:
        src = sources[name]
        sname, key = src["source_name"], src["config_key"]
        for step in steps_sel:
            sc = src["step_cfg"][step]
            w = str(sc["slots"])
            if step == "s1":
                argv = [python, str(SCRIPTS / "make_mcs_swath_masks.py"), "-c", str(CONFIG_DIR / src["mcs_config"]), "--workers", w]
                outs, deps = step_outputs(src, step, root, pcts), []
            elif step == "s2":
                if src["step2"] == "era5_imerg":
                    argv = [python, str(SCRIPTS / "combine_era5_imerg_tracking_masks.py"), "--era5_zarr", era5_zarr or defaults["era5_zarr"]]
                else:
                    argv = [python, str(SCRIPTS / "combine_tracking_masks.py"), "-c", cfg_sources, "--source", key]
                outs, deps = step_outputs(src, step, root, pcts), [(name, "s1")]
            elif step == "s3":
                argv = [python, str(SCRIPTS / "make_cooccurrence_masks.py"), "-c", cfg_sources, "--source", key, "--workers", w]
                outs, deps = step_outputs(src, step, root, pcts), [(name, "s2")]
            elif step == "monthly":
                argv = [python, str(SCRIPTS / "calc_monthly_rainmap_by_cof.py"), "-c", cfg_sources, "--source", key, "--workers", "12"]
                outs, deps = step_outputs(src, step, root, pcts), [(name, "s3")]
            elif step == "thresholds":
                in_zarr, in_var = _thresholds_input(src, root)
                argv = [python, str(SCRIPTS / "calc_extreme_precip_thresholds.py"), "--catalog_source", key, "--config_file", cfg_sources,
                        "--percentiles", *[str(p) for p in pcts], "--min_precip_threshold", str(defaults["min_precip_threshold"]),
                        "--cell_chunk_size", str(defaults["threshold_cell_chunk_size"]), "--output_dir", f"{root}extreme_precip/",
                        "--input_zarr", in_zarr]
                if in_var:
                    argv += ["--input_var", in_var]
                outs = step_outputs(src, step, root, pcts)
                deps = [(name, "s1")] if src.get("thresholds", {}).get("from_step1") else []
            else:  # attribution
                argv = [python, str(SCRIPTS / "calc_stormtype_extreme_precip_spatial.py"), "--catalog_source", key, "--config_file", cfg_sources,
                        "--percentiles", *[f"P{p}" for p in pcts], "--n_workers", w, "--output_dir", f"{root}extreme_precip"]
                outs = step_outputs(src, step, root, pcts)
                deps = [(name, "s3"), (name, "thresholds")]
            argv = argv + list((extra_args or {}).get(step, []))
            env = dict(env_base)
            env.update({k: str(v) for k, v in sc.get("env", {}).items()})
            tasks[(name, step)] = Task(name, step, argv, env, int(sc["slots"]), float(sc["est_min"]), float(sc["timeout_min"]), outs, deps)
    return tasks


# ---------------------------------------------------------------------------------------------------------------------------
# Markers and existing outputs
# ---------------------------------------------------------------------------------------------------------------------------
def marker_path(root, source, step):
    return Path(root) / "pipeline_state" / source / f"{step}.json"


def read_marker(root, source, step):
    try:
        with open(marker_path(root, source, step)) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def write_marker(root, task, status, git, extra=None):
    path = marker_path(root, task.source, task.step)
    path.parent.mkdir(parents=True, exist_ok=True)
    doc = {"source": task.source, "step": task.step, "status": status, "rc": task.rc,
           "start": _iso(task.start), "end": _iso(task.end), "minutes": round((task.end - task.start) / 60, 2) if task.start and task.end else None,
           "cmd": task.argv, "host": socket.gethostname(), "git": git, "data_root": root, "log": task.log_path,
           "peak_rss_gb": round(task.peak_rss_gb, 2), "outputs": task.outputs}
    doc.update(extra or {})
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump(doc, f, indent=1)
    os.replace(tmp, path)


def outputs_exist(task):
    return all(os.path.exists(o) for o in task.outputs)


def any_output_exists(task):
    return any(os.path.exists(o) for o in task.outputs)


def _iso(t):
    return datetime.datetime.fromtimestamp(t).isoformat(timespec="seconds") if t else None


def git_state():
    try:
        head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO, capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no", "--", "scripts", "src", "config"], cwd=REPO,
                               capture_output=True, text=True).stdout.strip()
        return head + (" (tracked files in scripts/src/config have uncommitted changes)" if dirty else "")
    except Exception:  # noqa: BLE001
        return "unknown"


# ---------------------------------------------------------------------------------------------------------------------------
# Plan: what runs, what is skipped, what conflicts
# ---------------------------------------------------------------------------------------------------------------------------
def topological_order(tasks):
    order, seen = [], set()

    def visit(k):
        if k in seen:
            return
        seen.add(k)
        for d in tasks[k].deps:
            if d in tasks:
                visit(d)
        order.append(k)
    # steps that are not in STEP_ORDER (another pipeline that reuses this engine, see run_etc_pipeline.py) sort after the known ones,
    # in the order the tasks were built; the dependencies alone decide when a step can start
    for k in sorted(tasks, key=lambda k: (STEP_ORDER.index(k[1]) if k[1] in STEP_ORDER else len(STEP_ORDER), k[0])):
        visit(k)
    return order


def prepare_plan(tasks, root, resume, force, outputs_of):
    """
    Decide which tasks run. Returns (conflicts, missing, warnings): tasks whose output exists without a marker of this runner
    (they would be overwritten), prerequisite steps that are neither selected nor present, and notes about prerequisites that
    exist without a marker. outputs_of((source, step)) gives the output paths of any step.
    """
    conflicts, missing, warnings = [], [], []
    for k in topological_order(tasks):
        t = tasks[k]
        marker = read_marker(root, t.source, t.step)
        upstream_reruns = any(d in tasks and tasks[d].will_run for d in t.deps)
        if resume and marker and marker.get("status") == "success" and outputs_exist(t) and not upstream_reruns:
            t.will_run, t.note, t.status = False, "done (marker)", "done"
        elif any_output_exists(t) and marker is None and not force:
            conflicts.append(t)
        for d in t.deps:
            if d in tasks:
                continue
            dm = read_marker(root, *d)
            if dm and dm.get("status") == "success":
                continue
            outs = outputs_of(d)
            if outs and all(os.path.exists(o) for o in outs):
                warnings.append(f"{t.label}: needs {d[0]}/{d[1]}, found its output without a marker of this runner (made outside it)")
                continue
            missing.append((t, d))
    return conflicts, missing, warnings


def critical_paths(tasks):
    """Remaining chain length in minutes from each task to the end of the run (its own duration included)."""
    succ = {k: [] for k in tasks}
    for k, t in tasks.items():
        for d in t.deps:
            if d in tasks:
                succ[d].append(k)
    memo = {}

    def cp(k):
        if k not in memo:
            memo[k] = tasks[k].est_min + max([cp(s) for s in succ[k]] or [0])
        return memo[k]
    return {k: cp(k) for k in tasks}


def simulate(tasks, max_slots, stagger_min):
    """Schedule with the estimated durations in virtual time; returns {key: (start, end)} in minutes and the makespan."""
    cp = critical_paths(tasks)
    todo = {k for k, t in tasks.items() if t.will_run}
    finish, times, running, used, now, last = {}, {}, [], 0, 0.0, -1e9
    for k, t in tasks.items():
        if not t.will_run:
            finish[k] = 0.0
    while todo or running:
        ready = sorted((k for k in todo if all(d not in tasks or d in finish for d in tasks[k].deps)), key=lambda k: -cp[k])
        launched = False
        for k in ready:
            t = tasks[k]
            if (used + t.slots <= max_slots or used == 0) and (now - last >= stagger_min - 1e-9 or used == 0):
                times[k] = (now, now + t.est_min)
                running.append((now + t.est_min, k))
                running.sort()
                used += t.slots
                todo.discard(k)
                last = now
                launched = True
                break
        if launched:
            continue
        if not running:
            break   # a task whose prerequisite is neither selected nor present cannot start; the plan check reports it
        nxt = running[0][0]
        if ready and now - last < stagger_min - 1e-9:   # same tolerance as the start test above
            nxt = min(nxt, last + stagger_min)
        if nxt <= now and not (running and running[0][0] <= now + 1e-9):
            raise RuntimeError("simulate: the virtual clock does not advance")   # cannot happen; better than a silent endless loop
        now = max(now, nxt)
        while running and running[0][0] <= now + 1e-9:
            end, k = running.pop(0)
            finish[k] = end
            used -= tasks[k].slots
    return times, max((e for _, e in times.values()), default=0.0)


# ---------------------------------------------------------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------------------------------------------------------
def preflight(tasks, sources, defaults, python, era5_zarr):
    """Checks before anything runs; returns (problems, notes). Reads directory listings and file names only."""
    problems, notes = [], []
    # the environment that runs the steps must have zarr 2 (zarr 3 writes another layout and breaks the readers)
    try:
        out = subprocess.run([python, "-c", "import zarr; print(zarr.__version__)"], capture_output=True, text=True, timeout=120)
        version = out.stdout.strip()
        if not version.startswith("2."):
            problems.append(f"{python}: zarr {version or '(not importable)'}; the pipeline needs zarr 2.x")
        else:
            notes.append(f"python {python}: zarr {version}")
    except Exception as exc:  # noqa: BLE001
        problems.append(f"cannot run {python}: {exc}")
    cfg_sources = yaml.safe_load(open(CONFIG_DIR / "config_sources.yaml"))
    stores = {}
    for (name, step), t in tasks.items():
        if not t.will_run:
            continue
        src = sources[name]
        if step == "s1":
            mcs_cfg = yaml.safe_load(open(CONFIG_DIR / src["mcs_config"]))
            base = mcs_cfg["zarr_output_presets"]["healpix"]["out_filebase"]
            stores[f"{mcs_cfg['root_path']}{mcs_cfg['pixel_path_name']}/{base}hp8_v1.zarr"] = f"{name}: hourly MCS masks (Step 1 input)"
            for extra in src.get("check_stores", []):
                stores[extra] = f"{name}: catalog input (Step 1)"
        if step == "s2":
            if src["step2"] == "era5_imerg":
                stores[era5_zarr or defaults["era5_zarr"]] = f"{name}: ERA5 AR/TC/ETC masks (Step 2 input)"
            else:
                c = cfg_sources[src["config_key"]]
                for kind in ("AR_tracks", "TC_test_tracks", "ETC_test_tracks"):
                    prefix = f"{kind}_{c['source_te']}_{c['source_res']}."
                    try:
                        n = sum(1 for f in os.listdir(c["dir_te"]) if f.startswith(prefix) and f.endswith(".nc"))
                    except OSError:
                        n = 0
                    (notes if n else problems).append(f"{name}: {n} {kind} files in {c['dir_te']}")
    for store, what in stores.items():
        try:
            scan = scan_store(store)
        except Exception as exc:  # noqa: BLE001
            problems.append(f"{what}: {store} cannot be checked ({type(exc).__name__}: {exc})")
            continue
        complete, lines = summarize(store, scan)
        if complete:
            notes.append(f"{what}: complete ({sum(r['expected'] for r in scan.values())} chunks)")
        else:
            problems.append(f"{what}: {store} is INCOMPLETE: " + "; ".join(t.strip() for p, t in lines if p))
    return problems, notes


# ---------------------------------------------------------------------------------------------------------------------------
# Resources
# ---------------------------------------------------------------------------------------------------------------------------
def mem_available_gb():
    if psutil:
        return psutil.virtual_memory().available / 1e9
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1e6
    return float("inf")


def mem_total_gb():
    if psutil:
        return psutil.virtual_memory().total / 1e9
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) / 1e6
    return 0.0


def tree_rss_gb(pid):
    """Resident memory of a process and all its children (Dask workers), and how many processes that is."""
    if not psutil:
        return 0.0, 0
    try:
        procs = [psutil.Process(pid)] + psutil.Process(pid).children(recursive=True)
    except psutil.Error:
        return 0.0, 0
    total = 0
    for p in procs:
        try:
            total += p.memory_info().rss
        except psutil.Error:
            pass
    return total / 1e9, len(procs)


# ---------------------------------------------------------------------------------------------------------------------------
# The scheduler
# ---------------------------------------------------------------------------------------------------------------------------
class Runner:
    def __init__(self, tasks, root, log_dir, max_slots, stagger_sec, min_free_gb, git, poll_sec=1.0, sample_sec=30.0):
        self.tasks, self.root, self.log_dir, self.git = tasks, root, Path(log_dir), git
        self.max_slots, self.stagger, self.min_free_gb = max_slots, stagger_sec, min_free_gb
        self.poll, self.sample_every = poll_sec, sample_sec
        self.cp = critical_paths(tasks)
        self.used = 0
        self.max_used = 0
        self.min_avail = float("inf")
        self.stop_signal = None
        self.t0 = time.time()
        self.last_launch = 0.0
        self.last_sample = 0.0
        self.last_status = 0.0
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.resources = open(self.log_dir / "resources.csv", "w")
        self.resources.write("time,elapsed_min,task,rss_gb,processes,node_mem_available_gb,slots_used,running\n")

    # -- helpers ------------------------------------------------------------------------------------------------------------
    def say(self, text):
        print(f"{datetime.datetime.now():%H:%M:%S} [{(time.time() - self.t0) / 60:6.1f} min] {text}", flush=True)

    def ready(self, t):
        return all(d not in self.tasks or self.tasks[d].status in ("ok", "done") for d in t.deps)

    def blocked_by_failure(self, t):
        return any(d in self.tasks and self.tasks[d].status in ("failed", "skipped", "interrupted") for d in t.deps)

    def running(self):
        return [t for t in self.tasks.values() if t.status == "running"]

    def write_status(self):
        st = {"updated": _iso(time.time()), "elapsed_min": round((time.time() - self.t0) / 60, 1), "slots_used": self.used,
              "max_slots": self.max_slots, "mem_available_gb": round(mem_available_gb(), 1),
              "running": [t.label for t in self.running()],
              "done": [t.label for t in self.tasks.values() if t.status in ("ok", "done")],
              "failed": [t.label for t in self.tasks.values() if t.status in ("failed", "interrupted")],
              "skipped": [t.label for t in self.tasks.values() if t.status == "skipped"],
              "pending": [t.label for t in self.tasks.values() if t.status == "pending"]}
        with open(self.log_dir / "status.json", "w") as f:
            json.dump(st, f, indent=1)

    # -- launching and finishing -------------------------------------------------------------------------------------------
    def launch(self, t):
        t.log_path = str(self.log_dir / f"{t.source}_{t.step}.log")
        env = dict(os.environ)
        env["COF_DATA_ROOT"] = self.root             # always the runner's root, whatever the caller's environment says
        env.update(t.env)                            # a task may set its own (the ETC pipeline reads the COF products from another root)
        logf = open(t.log_path, "wb")
        t.start = time.time()
        t.proc = subprocess.Popen(t.argv, cwd=REPO, env=env, stdout=logf, stderr=subprocess.STDOUT, start_new_session=True)
        logf.close()
        t.status = "running"
        self.used += t.slots
        self.max_used = max(self.max_used, self.used)
        self.last_launch = time.time()
        write_marker(self.root, t, "running", self.git)
        self.say(f"start   {t.label:34s} slots {t.slots:3d} (in use {self.used}/{self.max_slots}), estimated {t.est_min:.0f} min, log {t.log_path}")

    def finish(self, t, rc, status=None):
        t.end, t.rc = time.time(), rc
        t.status = status or ("ok" if rc == 0 else "failed")
        self.used -= t.slots
        write_marker(self.root, t, {"ok": "success", "failed": "failed", "interrupted": "interrupted"}[t.status], self.git)
        self.say(f"{'done ' if t.status == 'ok' else t.status.upper():7s} {t.label:34s} rc={rc} after {(t.end - t.start) / 60:.1f} min, peak {t.peak_rss_gb:.1f} GB")

    def terminate(self, t, grace=20):
        """Stop a task and all its processes (its own process group)."""
        try:
            os.killpg(t.proc.pid, signal.SIGTERM)
        except (ProcessLookupError, PermissionError):
            pass
        deadline = time.time() + grace
        while t.proc.poll() is None and time.time() < deadline:
            time.sleep(0.2)
        if t.proc.poll() is None:
            try:
                os.killpg(t.proc.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
            t.proc.wait()

    def sample(self):
        avail = mem_available_gb()
        self.min_avail = min(self.min_avail, avail)
        running = self.running()
        now = f"{datetime.datetime.now():%H:%M:%S}"
        el = f"{(time.time() - self.t0) / 60:.2f}"
        for t in running:
            rss, n = tree_rss_gb(t.proc.pid)
            t.peak_rss_gb = max(t.peak_rss_gb, rss)
            self.resources.write(f"{now},{el},{t.label},{rss:.2f},{n},{avail:.1f},{self.used},{len(running)}\n")
        self.resources.write(f"{now},{el},NODE,,,{avail:.1f},{self.used},{len(running)}\n")
        self.resources.flush()

    # -- the loop ----------------------------------------------------------------------------------------------------------
    def run(self):
        for sig in (signal.SIGTERM, signal.SIGINT):
            signal.signal(sig, lambda s, f: setattr(self, "stop_signal", s))
        for t in self.tasks.values():
            if not t.will_run:
                t.status = "done"
        self.say(f"{sum(t.will_run for t in self.tasks.values())} steps to run, {sum(not t.will_run for t in self.tasks.values())} already done; "
                 f"slot budget {self.max_slots}, node memory {mem_total_gb():.0f} GB")
        while True:
            # finished steps
            for t in self.running():
                rc = t.proc.poll()
                if rc is None and time.time() - t.start > t.timeout_min * 60:
                    self.say(f"TIMEOUT {t.label} after {t.timeout_min:.0f} min: stopping it")
                    self.terminate(t)
                    rc = t.proc.returncode if t.proc.returncode not in (None, 0) else -1
                if rc is not None:
                    self.finish(t, rc)
            if self.stop_signal:
                break
            # steps that cannot run any more
            for t in self.tasks.values():
                if t.status == "pending" and t.will_run and self.blocked_by_failure(t):
                    t.status = "skipped"
                    t.note = "a step it needs did not finish"
                    self.say(f"skipped {t.label} (a step it needs did not finish)")
            pending = [t for t in self.tasks.values() if t.status == "pending" and t.will_run]
            if not pending and not self.running():
                break
            # new steps: longest remaining chain first, within the slot budget, staggered, only with memory to spare
            ready = sorted((t for t in pending if self.ready(t)), key=lambda t: -self.cp[t.key])
            for t in ready:
                fits = self.used + t.slots <= self.max_slots or self.used == 0
                spaced = time.time() - self.last_launch >= self.stagger or self.used == 0
                if fits and spaced:
                    if self.running() and mem_available_gb() < self.min_free_gb:
                        break   # wait for memory; something is running, so it will free up
                    self.launch(t)
            now = time.time()
            if now - self.last_sample >= self.sample_every:
                self.sample()
                self.last_sample = now
            if now - self.last_status >= 5:
                self.write_status()
                self.last_status = now
            if now - getattr(self, "_last_report", 0) >= 300:
                self._last_report = now
                self.say("running: " + (", ".join(f"{t.label} ({(now - t.start) / 60:.0f} min)" for t in self.running()) or "none") +
                         f" | slots {self.used}/{self.max_slots} | node memory available {mem_available_gb():.0f} GB")
            time.sleep(self.poll)
        if self.stop_signal:
            self.say(f"signal {self.stop_signal} received: stopping {len(self.running())} running step(s); the run can be continued with --resume")
            for t in self.running():
                self.terminate(t)
                self.finish(t, t.proc.returncode, "interrupted")
        self.sample()
        self.write_status()
        self.resources.close()
        return self.summary()

    # -- report -----------------------------------------------------------------------------------------------------------
    def summary(self):
        rows = []
        for k in topological_order(self.tasks):
            t = self.tasks[k]
            rows.append((t.source, t.step, t.status if t.will_run else "done (marker)",
                         f"{(t.start - self.t0) / 60:6.1f}" if t.start else "     -",
                         f"{(t.end - t.start) / 60:6.1f}" if t.end and t.start else "     -", t.slots, f"{t.peak_rss_gb:6.1f}" if t.start else "     -"))
        lines = [f"{'source':24s} {'step':12s} {'status':14s} {'start min':>9s} {'minutes':>8s} {'slots':>5s} {'peak RSS GB':>11s}"]
        lines += [f"{r[0]:24s} {r[1]:12s} {r[2]:14s} {r[3]:>9s} {r[4]:>8s} {r[5]:5d} {r[6]:>11s}" for r in rows]
        n = {s: sum(1 for t in self.tasks.values() if (t.status if t.will_run else "done") == s) for s in ("ok", "done", "failed", "skipped", "interrupted", "pending")}
        lines.append("")
        lines.append(f"steps: {n['ok']} ok, {n['done']} already done, {n['failed']} failed, {n['skipped']} skipped after a failure, "
                     f"{n['interrupted']} interrupted, {n['pending']} not started; wall time {(time.time() - self.t0) / 60:.1f} min; "
                     f"most slots in use {self.max_used}/{self.max_slots}; lowest node memory available {self.min_avail:.0f} GB")
        good = n["failed"] == 0 and n["skipped"] == 0 and n["interrupted"] == 0 and n["pending"] == 0
        return good, "\n".join(lines)


# ---------------------------------------------------------------------------------------------------------------------------
# Self-test on fake steps
# ---------------------------------------------------------------------------------------------------------------------------
def _fake_tasks(spec, out_dir):
    """Tasks that only sleep and touch their output: [{'source','step','seconds','rc','slots','deps':[[source, step], ...]}]"""
    code = "import sys,time,os; time.sleep(float(sys.argv[2])); os.makedirs(os.path.dirname(sys.argv[1]), exist_ok=True); open(sys.argv[1],'w').write('x'); sys.exit(int(sys.argv[3]))"
    tasks = {}
    for s in spec:
        out = os.path.join(out_dir, f"{s['source']}_{s['step']}.out")
        tasks[(s["source"], s["step"])] = Task(s["source"], s["step"], [sys.executable, "-c", code, out, str(s["seconds"]), str(s.get("rc", 0))], {},
                                                 s.get("slots", 1), s["seconds"] / 60, 5.0, [out], [tuple(d) for d in s.get("deps", [])])
    return tasks


def self_test():
    """Order, slot budget, failure skipping, resume, overwrite protection and interruption, on fake steps that only sleep."""
    fails = []

    def check(name, ok, detail=""):
        print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
        if not ok:
            fails.append(name)

    tmp = tempfile.mkdtemp(prefix="cof_selftest_")
    root, out = f"{tmp}/root/", f"{tmp}/out"
    outputs_of = lambda d: [os.path.join(out, f"{d[0]}_{d[1]}.out")]   # noqa: E731
    spec = []
    for src, s1_sec in (("A", 1.0), ("B", 0.6), ("C", 0.8)):
        spec += [dict(source=src, step="s1", seconds=s1_sec, slots=4), dict(source=src, step="s2", seconds=0.3, slots=2, deps=[[src, "s1"]]),
                 dict(source=src, step="s3", seconds=0.5, slots=4, deps=[[src, "s2"]]), dict(source=src, step="monthly", seconds=0.3, slots=3, deps=[[src, "s3"]]),
                 dict(source=src, step="thresholds", seconds=0.4, slots=2, deps=[[src, "s1"]]),
                 dict(source=src, step="attribution", seconds=0.4, slots=3, deps=[[src, "s3"], [src, "thresholds"]])]
    # 1. a normal run: order and slot budget
    tasks = _fake_tasks(spec, out)
    prepare_plan(tasks, root, resume=False, force=False, outputs_of=outputs_of)
    runner = Runner(tasks, root, f"{tmp}/log1", max_slots=9, stagger_sec=0.05, min_free_gb=0, git="test", poll_sec=0.05, sample_sec=0.2)
    good, text = runner.run()
    check("normal run: everything ok", good and all(t.status == "ok" for t in tasks.values()))
    order_ok = all(tasks[d].end <= t.start + 1e-3 for t in tasks.values() for d in t.deps)
    check("normal run: every step started after the steps it needs had finished", order_ok)
    check("normal run: slots in use never exceeded the budget", runner.max_used <= 9, f"max {runner.max_used} of 9")
    check("normal run: steps of different sources overlapped", any(a.start < b.end and b.start < a.end for a in tasks.values() for b in tasks.values() if a.source != b.source))
    check("normal run: markers written", all(read_marker(root, t.source, t.step)["status"] == "success" for t in tasks.values()))
    # 2. overwrite protection
    tasks = _fake_tasks(spec, out)
    for p in Path(root).glob("pipeline_state/*/*.json"):
        p.unlink()
    conflicts, missing, _ = prepare_plan(tasks, root, resume=False, force=False, outputs_of=outputs_of)
    check("outputs without a marker are reported as conflicts (not overwritten without --force)", len(conflicts) == len(spec), f"{len(conflicts)} conflicts")
    conflicts, _, _ = prepare_plan(_fake_tasks(spec, out), root, resume=False, force=True, outputs_of=outputs_of)
    check("--force accepts them", not conflicts)
    # 3. a failure skips only its own dependents; resume finishes the run
    shutil_out = Path(out)
    for f in shutil_out.glob("*"):
        f.unlink()
    shutil_root = Path(root)
    for p in shutil_root.glob("pipeline_state/*/*.json"):
        p.unlink()
    spec_fail = copy.deepcopy(spec)
    for s in spec_fail:
        if s["source"] == "B" and s["step"] == "s2":
            s["rc"] = 3
    tasks = _fake_tasks(spec_fail, out)
    prepare_plan(tasks, root, resume=False, force=False, outputs_of=outputs_of)
    runner = Runner(tasks, root, f"{tmp}/log2", max_slots=9, stagger_sec=0.05, min_free_gb=0, git="test", poll_sec=0.05, sample_sec=0.2)
    good, text = runner.run()
    check("failure: the run reports it", not good)
    check("failure: the failed step is failed, its dependents (s3, monthly, attribution) are skipped",
          tasks[("B", "s2")].status == "failed" and all(tasks[("B", s)].status == "skipped" for s in ("s3", "monthly", "attribution")))
    check("failure: the other sources and the independent step (thresholds) finished",
          all(t.status == "ok" for t in tasks.values() if t.source != "B") and tasks[("B", "thresholds")].status == "ok")
    tasks = _fake_tasks(spec, out)   # the step works now
    conflicts, missing, _ = prepare_plan(tasks, root, resume=True, force=False, outputs_of=outputs_of)
    check("resume: no conflicts, nothing missing", not conflicts and not missing)
    check("resume: steps with a success marker are skipped, the failed step and what follows are planned",
          all(not tasks[k].will_run for k in tasks if k[0] != "B") and not tasks[("B", "s1")].will_run and not tasks[("B", "thresholds")].will_run
          and all(tasks[("B", s)].will_run for s in ("s2", "s3", "monthly", "attribution")))
    runner = Runner(tasks, root, f"{tmp}/log3", max_slots=9, stagger_sec=0.05, min_free_gb=0, git="test", poll_sec=0.05, sample_sec=0.2)
    good, text = runner.run()
    check("resume: completes", good and all(t.status in ("ok", "done") for t in tasks.values()))
    # 4. prerequisites that are not selected
    for f in Path(out).glob("A_s3.out"):   # the selected steps' own outputs must not exist for this scenario
        f.unlink()
    for f in Path(out).glob("A_monthly.out"):
        f.unlink()
    tasks = _fake_tasks([s for s in spec if s["source"] == "A" and s["step"] in ("s3", "monthly")], out)
    conflicts, missing, warns = prepare_plan(tasks, f"{tmp}/other_root/", resume=False, force=False, outputs_of=lambda d: [f"{tmp}/none/{d[0]}_{d[1]}.out"])
    check("missing prerequisites are reported (s3 needs s2, which is neither selected nor present)", any(d == ("A", "s2") for _, d in missing), str([(t.label, d) for t, d in missing]))
    tasks = _fake_tasks([s for s in spec if s["source"] == "A" and s["step"] in ("s3", "monthly")], out)
    conflicts, missing, warns = prepare_plan(tasks, f"{tmp}/other_root/", resume=False, force=False, outputs_of=outputs_of)
    check("a prerequisite that exists without a marker is accepted with a warning", not missing and len(warns) == 1, str(warns))
    print(f"\nRESULT: {'PASS' if not fails else 'FAIL ' + str(fails)}")
    return 0 if not fails else 1


def interrupt_test():
    """SIGTERM to a running runner: children stop, the markers say 'interrupted', the exit status is 1."""
    tmp = tempfile.mkdtemp(prefix="cof_selftest_int_")
    spec = [dict(source="A", step="s1", seconds=30, slots=2), dict(source="A", step="s2", seconds=1, slots=1, deps=[["A", "s1"]])]
    spec_file = f"{tmp}/spec.json"
    json.dump(spec, open(spec_file, "w"))
    proc = subprocess.Popen([sys.executable, __file__, "--fake-spec", spec_file, "--data-root", f"{tmp}/root", "--stagger-sec", "0", "--max-slots", "4"],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    deadline = time.time() + 20
    while time.time() < deadline and not Path(f"{tmp}/root/pipeline_state/A/s1.json").exists():
        time.sleep(0.2)
    time.sleep(1.0)
    proc.send_signal(signal.SIGTERM)
    try:
        proc.wait(timeout=40)
    except subprocess.TimeoutExpired:
        proc.kill()
    m = read_marker(f"{tmp}/root", "A", "s1")
    ok = proc.returncode == 1 and m and m["status"] == "interrupted" and read_marker(f"{tmp}/root", "A", "s2") is None
    print(f"[{'PASS' if ok else 'FAIL'}] interrupt: exit status {proc.returncode}, marker of the running step '{m and m['status']}', the dependent step never started")
    return 0 if ok else 1


# ---------------------------------------------------------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Run Analyses 1 and 2 of the COF pipeline for any set of sources, in dependency order",
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__.split("Examples")[1].split("Exit status")[0])
    p.add_argument("--data-root", required=False, help="root of all products (mcs_masks/, all_masks/, cof_masks/, extreme_precip/); use a test area, "
                   "production is /pscratch/sd/w/wcmca1/hackathon/")
    p.add_argument("--sources", nargs="*", default=None, help="sources to run, names or aliases (default: all in the registry)")
    p.add_argument("--analysis", choices=["1", "2", "both"], default="both", help="1: total precipitation (s1 s2 s3 monthly), 2: extreme precipitation "
                   "(s1 s2 s3 thresholds attribution), both (default)")
    p.add_argument("--steps", nargs="*", default=None, help=f"only these steps ({', '.join(STEP_ORDER)}; 1 2 3 thr attr accepted); the steps they need must exist")
    p.add_argument("--from", dest="from_step", default=None, help="only this step and the ones after it in the order s1 s2 s3 monthly thresholds attribution")
    p.add_argument("--registry", default=str(CONFIG_DIR / "config_pipeline.yaml"))
    p.add_argument("--python", default=None, help="python that runs the steps (default: the registry's, the hackathon environment)")
    p.add_argument("--max-slots", type=int, default=None, help="CPU slots (workers x threads) that may be busy at once (default: 80%% of the logical CPUs)")
    p.add_argument("--stagger-sec", type=float, default=20.0, help="minimum seconds between two step starts (default 20)")
    p.add_argument("--min-free-gb", type=float, default=60.0, help="do not start a step while less memory than this is available (default 60 GB)")
    p.add_argument("--resume", action="store_true", help="skip steps whose marker says success; re-run the rest (also a step that was interrupted)")
    p.add_argument("--force", action="store_true", help="overwrite outputs that this runner did not create")
    p.add_argument("--dry-run", action="store_true", help="show the graph, the commands and the estimated schedule; run nothing")
    p.add_argument("--preflight-only", action="store_true", help="run the input checks and stop")
    p.add_argument("--skip-preflight", action="store_true", help="do not check the inputs first")
    p.add_argument("--step-args", nargs=2, action="append", metavar=("STEP", "ARGS"), default=[],
                   help='extra arguments appended to the command of one step, e.g. --step-args s1 "--test-steps 24 --workers 4" (repeatable); '
                        'for tests, or a special run')
    p.add_argument("--era5-zarr", default=None, help="ERA5 AR/TC/ETC store for IMERG's Step 2 (default: the registry's, the CFS copy)")
    p.add_argument("--log-dir", default=None, help="logs, status.json and resources.csv (default: <data-root>/pipeline_logs/<run id>)")
    p.add_argument("--self-test", action="store_true", help="test the scheduler on fake steps and exit")
    p.add_argument("--fake-spec", default=None, help=argparse.SUPPRESS)   # used by the interrupt test
    return p.parse_args()


def main():
    args = parse_args()
    if args.self_test:
        rc = self_test()
        rc |= interrupt_test()
        sys.exit(rc)
    if not args.data_root:
        sys.exit("--data-root is required (a test area, or /pscratch/sd/w/wcmca1/hackathon/ for production)")
    root = args.data_root.rstrip("/") + "/"
    logical = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    max_slots = args.max_slots or max(1, int(0.8 * logical))
    run_id = time.strftime("%Y%m%d_%H%M%S")
    log_dir = args.log_dir or f"{root}pipeline_logs/{run_id}"

    if args.fake_spec:   # the interrupt test: fake steps, the real scheduler
        spec = json.load(open(args.fake_spec))
        tasks = _fake_tasks(spec, os.path.join(os.path.dirname(args.fake_spec), "out"))
        prepare_plan(tasks, root, resume=False, force=False, outputs_of=lambda d: [])
        good, text = Runner(tasks, root, log_dir, max_slots, args.stagger_sec, 0, "test", poll_sec=0.1, sample_sec=1).run()
        print(text)
        sys.exit(0 if good else 1)

    defaults, sources = load_registry(args.registry)
    python = args.python or defaults["python"]
    names = resolve_sources(args.sources, sources)
    steps = normalize_steps(args.steps) if args.steps else list(ANALYSIS_STEPS[args.analysis])
    if args.from_step:
        first = normalize_steps([args.from_step])[0]
        steps = [s for s in steps if STEP_ORDER.index(s) >= STEP_ORDER.index(first)]
    steps = [s for s in STEP_ORDER if s in steps]
    extra = {}
    for step_name, step_args in args.step_args:
        extra.setdefault(normalize_steps([step_name])[0], []).extend(shlex.split(step_args))
    tasks = build_tasks(names, steps, sources, defaults, root, python, args.era5_zarr, extra)
    conflicts, missing, warnings = prepare_plan(tasks, root, args.resume, args.force,
                                                lambda d: step_outputs(sources[d[0]], d[1], root, defaults["percentiles"]))

    print(f"Data root: {root}\nSources: {', '.join(names)}\nSteps: {', '.join(steps)}\nPython: {python}\nSlot budget: {max_slots} of {logical} logical CPUs")
    for w in warnings:
        print(f"WARNING: {w}")
    if missing:
        print("\nA step needs an earlier step that is neither selected nor present:")
        for t, d in missing:
            print(f"  {t.label} needs {d[0]}/{d[1]}")
        print("Add the steps to the selection (--steps / --analysis / --from) or run them first.")
        sys.exit(2)
    if conflicts:
        print("\nThese outputs exist and were not made by this runner (no marker); nothing is overwritten unless --force is given:")
        for t in conflicts:
            print(f"  {t.label}: {[o for o in t.outputs if os.path.exists(o)][0]}")
        sys.exit(2)

    if args.dry_run:
        times, span = simulate(tasks, max_slots, args.stagger_sec / 60)
        print("\nGraph (start and end are estimates in minutes from the start; steps run in parallel where their needs allow):")
        for k in topological_order(tasks):
            t = tasks[k]
            when = f"{times[k][0]:6.1f} -> {times[k][1]:6.1f}" if k in times else "   already done"
            print(f"  {t.label:34s} slots {t.slots:3d}  {when}   needs: {', '.join(f'{d[0]}/{d[1]}' for d in t.deps) or '-'}")
            print(f"      {' '.join(t.argv)}")
        print(f"\nEstimated wall time: {span:.0f} min ({span / 60:.2f} h) if the estimates hold and the node is not slowed down by the load.")
        return

    if not args.skip_preflight:
        problems, notes = preflight(tasks, sources, defaults, python, args.era5_zarr)
        for n in notes:
            print(f"  ok: {n}")
        for pr in problems:
            print(f"  PROBLEM: {pr}")
        if problems:
            print("Preflight failed; fix the inputs or use --skip-preflight.")
            sys.exit(2)
        print("Preflight passed.")
    if args.preflight_only:
        return

    git = git_state()
    print(f"Code version: {git}\nLogs: {log_dir}  (status.json is updated every few seconds)\n", flush=True)
    good, text = Runner(tasks, root, log_dir, max_slots, args.stagger_sec, args.min_free_gb, git).run()
    print("\n" + text)
    (Path(log_dir) / "summary.txt").write_text(text + "\n")
    sys.exit(0 if good else 1)


if __name__ == "__main__":
    main()
