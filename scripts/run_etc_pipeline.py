#!/usr/bin/env python
"""
Run Analysis 3 (the ETC composites) for any set of sources, in dependency order, on one node.

Per source (see config/config_etc_pipeline.yaml and docs/procedures/run_etc_pipeline.md):

    etc_cof                                            ETC tracks + the Step 3 overlap flags        -> <root>/etc_tracks/<src>_etc_cof_data.parquet
    pr                                                 precipitation around every ETC point         -> <root>/etc_data/<src>/single_vars/etc_2d_pr_all_all.zarr
    mask_<variable>  (seven COF mask variables)        COF masks around every ETC point             -> ... etc_2d_<variable>_all_all.zarr
    link_env                                           the unchanged environment stores, linked from etc_data/ (no copy, checked first)
    combine     (needs the four above)                 one multi-variable store                     -> <root>/etc_data/<src>/etc_2d_combined_all_all.zarr
    composites, stats   (need combine)                 composites and spatial statistics            -> <root>/etc_data/stats/

It reuses the engine of run_cof_pipeline.py (CPU-slot budget, memory gate, markers, --resume, protection of outputs that it did not create,
logs and status.json). The COF products (cof_masks/) and the track files are inputs, read from --cof-root (default: the production tree) and
the registry, not from --data-root, which receives the new products and is required: use a test area first.

Examples:
  python scripts/run_etc_pipeline.py --data-root /pscratch/sd/w/wcmca1/hackathon/tmp/etc_round1 --dry-run
  python scripts/run_etc_pipeline.py --data-root DIR --sources um scream          # names or aliases from the registry
  python scripts/run_etc_pipeline.py --data-root DIR --steps combine composites stats --resume
  python scripts/run_etc_pipeline.py --data-root DIR --step-args pr "--start_date 2020-03-01 --end_date 2020-03-08"    # a short test

Exit status: 0 all steps done; 1 a step failed, was skipped after a failure, or the run was interrupted; 2 usage error, refused
overwrite, missing prerequisite or failed preflight.

Author: Zhe Feng | zhe.feng@pnnl.gov
"""
import argparse
import contextlib
import io
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"
EXTRACT_DIR = REPO / "extract_environments"
CONFIG_DIR = REPO / "config"
for p in (SCRIPTS, EXTRACT_DIR, REPO):
    sys.path.insert(0, str(p))
import run_cof_pipeline as eng                          # noqa: E402  the engine (Task, Runner, prepare_plan, simulate, markers)
import submit_etc_extraction_jobs as sub                # noqa: E402  track files, catalog settings, COF mask variables, expected gaps
from check_zarr_store import scan_store, summarize      # noqa: E402

STEP_NAMES = ["etc_cof", "pr", "mask", "link_env", "combine", "composites", "stats"]
MASK_VARS = list(sub.COF_MASK_VARS)
COMPOSITE_FILES = [f"etc_2d_composite_{h}_{o}.nc" for h in ("nh", "sh") for o in ("all", "isolated", "mcs_only", "ar_only", "3way")]


# ---------------------------------------------------------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------------------------------------------------------
def load_registry(path):
    """Returns (defaults, {source: dict with 'step_cfg': {step: {slots, est_min, timeout_min}}})."""
    with open(path) as f:
        reg = yaml.safe_load(f)
    defaults, steps, sources = reg["defaults"], reg["steps"], {}
    for name, src in reg["sources"].items():
        src = dict(src, name=name)
        cfg = {}
        for step in STEP_NAMES:
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
        for alias in [name] + list(src.get("aliases", [])) + [src["submit_key"]]:
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
        if n.lower() not in STEP_NAMES:
            raise SystemExit(f"Unknown step '{n}'. Steps: {', '.join(STEP_NAMES)}")
        out.append(n.lower())
    return out


# ---------------------------------------------------------------------------------------------------------------------------
# Tasks
# ---------------------------------------------------------------------------------------------------------------------------
def source_info(src):
    """(model config of submit_etc_extraction_jobs, the source name of the scripts, the name of its COF store)."""
    mc = sub.MODEL_CONFIGS[src["submit_key"]]
    a3 = Path(mc["output_dir"]).parent.name              # scream, icon_d3hp003, ..., era5
    return mc, a3, ("IMERGv7" if a3 == "era5" else a3)


def find_group(mc, kind):
    """The extraction group of a model config: 'pr' (the one that extracts pr) or 'cof_mask'."""
    for g in mc["job_groups"]:
        if kind == "cof_mask" and g.get("cof_mask"):
            return g
        if kind == "pr" and not g.get("cof_mask") and "pr" in g["variables"]:
            return g
    raise SystemExit(f"no '{kind}' group in the extraction settings of {mc['job_name']}")


def extract_argv(python, mc, group, variable, out_dir):
    """The command of extract_etc_2d_vars.py for one variable, built by the same code that makes the Slurm scripts, with another output folder."""
    # the text is made for a shell script (backslash-newline continuations); shlex would keep each escaped newline as a token of its own
    toks = shlex.split(sub.build_python_cmd(mc, group).replace("\\\n", " "))
    toks[0] = python
    toks = [variable if t == "$CURRENT_VAR" else t for t in toks]
    toks[toks.index("--output_dir") + 1] = out_dir
    return toks


def expected_points(mc):
    """Number of storm points the extraction keeps: valid positions within 90 - radius degrees of the equator (as extract_etc_2d_vars.py)."""
    from src.env_extract_utilities import parse_etc_track_file
    with contextlib.redirect_stdout(io.StringIO()):
        df = parse_etc_track_file(mc["track_file"], unstructured_mesh=not mc.get("structured_mesh", False))
    radius = float(sub.RADIUS)
    return int(df.dropna(subset=["lat", "lon"])["lat"].between(-90 + radius, 90 - radius).sum())


def build_tasks(names, sources, defaults, root, python, cof_root, env_from, extra):
    """The full graph for the given sources (all steps): {(source, step): Task}, and per-source information."""
    tasks, infos = {}, {}
    base_env = dict(defaults["env"])
    for name in names:
        src = sources[name]
        mc, a3, cof = source_info(src)
        cfgs = src["step_cfg"]
        single = f"{root}etc_data/{a3}/single_vars/"
        n_pts = expected_points(mc)
        env_exclude = ["pr", *MASK_VARS, *src.get("env_exclude", [])]     # not linked: extracted again, or stale (see the registry)
        infos[name] = {"mc": mc, "a3": a3, "cof": cof, "single": single, "expected_points": n_pts, "env_exclude": env_exclude}

        def add(step, cfg_key, argv, outputs, deps, env=None):
            c = cfgs[cfg_key]
            tasks[(name, step)] = eng.Task(name, step, argv + extra.get(cfg_key, []), {**base_env, **(env or {})}, int(c["slots"]), float(c["est_min"]),
                                           float(c["timeout_min"]), outputs, [(name, d) for d in deps])

        add("etc_cof", "etc_cof",
            [python, str(SCRIPTS / "combine_etc_cof_data.py"), "--source", a3, "--etc_dir", str(Path(mc["track_file"]).parent) + "/",
             "--cof_dir", f"{cof_root}cof_masks/stats/", "--output_dir", f"{root}etc_tracks/"],
            [f"{root}etc_tracks/{a3}_etc_cof_data.parquet"], [])
        cof_env = {"COF_DATA_ROOT": cof_root}
        add("pr", "pr", extract_argv(python, mc, find_group(mc, "pr"), "pr", single), [f"{single}etc_2d_pr_all_all.zarr"], [], cof_env)
        mask_group = find_group(mc, "cof_mask")
        for v in MASK_VARS:
            add(f"mask_{v}", "mask", extract_argv(python, mc, mask_group, v, single), [f"{single}etc_2d_{v}_all_all.zarr"], [], cof_env)
        add("link_env", "link_env",
            [python, str(EXTRACT_DIR / "link_etc_env_stores.py"), "--src-dir", f"{env_from}{a3}/single_vars", "--dst-dir", single,
             "--exclude", *env_exclude, "--expected-points", str(n_pts), "--min-stores", str(defaults["min_env_stores"])], [], [])
        add("combine", "combine",
            [python, str(EXTRACT_DIR / "combine_etc_2d_vars.py"), "--source", a3, "--input_dir", single, "--output_dir", f"{root}etc_data/{a3}",
             "--etc-path", f"{root}etc_tracks/"],
            [f"{root}etc_data/{a3}/etc_2d_combined_all_all.zarr"], ["etc_cof", "pr", *[f"mask_{v}" for v in MASK_VARS], "link_env"])
        add("composites", "composites",
            [python, str(SCRIPTS / "create_etc_composites.py"), "--source", a3, "--zarr-path", f"{root}etc_data/", "--out-dir", f"{root}etc_data/stats/{a3}"],
            [f"{root}etc_data/stats/{a3}/{f}" for f in COMPOSITE_FILES], ["combine"])
        add("stats", "stats",
            [python, str(SCRIPTS / "calc_etc_spatial_stats.py"), "--source", a3, "--zarr-path", f"{root}etc_data", "--output-dir", f"{root}etc_data/stats",
             "--basic-radius", "10.0", "--mask-radii", "10.0", "15.0", "--pr-radius", "10.0"],
            [f"{root}etc_data/stats/etc_spatial_stats_{a3}.nc"], ["combine"])
    return tasks, infos


def select(all_tasks, steps):
    """The tasks of the chosen steps ('mask' stands for the seven mask_<variable> tasks), in build order."""
    keep = set(steps) - {"mask"} | ({f"mask_{v}" for v in MASK_VARS} if "mask" in steps else set())
    return {k: t for k, t in all_tasks.items() if k[1] in keep}


# ---------------------------------------------------------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------------------------------------------------------
def preflight(tasks, infos, defaults, python, cof_root, env_from):
    """Checks before anything runs; returns (problems, notes). Reads directory listings and file names (chunk-key scans), no data."""
    problems, notes = [], []
    try:
        out = subprocess.run([python, "-c", "import zarr; print(zarr.__version__)"], capture_output=True, text=True, timeout=120)
        version = out.stdout.strip()
        (notes if version.startswith("2.") else problems).append(f"python {python}: zarr {version or '(not importable)'}" + ("" if version.startswith("2.") else "; the pipeline needs zarr 2.x"))
    except Exception as exc:  # noqa: BLE001
        problems.append(f"cannot run {python}: {exc}")
    seen = {}

    def check(store, what, points=None):
        if store not in seen:
            try:
                scan = scan_store(store)
                complete, lines = summarize(store, scan)
                seen[store] = "" if complete else "is INCOMPLETE: " + "; ".join(t.strip() for p_, t in lines if p_)
            except Exception as exc:  # noqa: BLE001
                seen[store] = f"cannot be checked ({type(exc).__name__}: {exc})"
        if seen[store]:
            problems.append(f"{what}: {store} {seen[store]}")
            return False
        if points is not None:
            import zarr
            n = int(zarr.open(store, mode="r")["time"].shape[0]) if "time" in zarr.open(store, mode="r") else None
            if n is not None and n != points:
                problems.append(f"{what}: {store} has {n} storm points, the current track file gives {points}")
                return False
        return True

    for name, info in infos.items():
        steps = {t.step for (n, _), t in tasks.items() if n == name and t.will_run}
        if not steps:
            continue
        mc, a3, cof, n_pts = info["mc"], info["a3"], info["cof"], info["expected_points"]
        if "etc_cof" in steps:
            for path, what in ((mc["track_file"], "ETC track file"), (f"{cof_root}cof_masks/stats/{cof}_etc_overlap_tracking.parquet", "Step 3 overlap parquet")):
                (notes if os.path.exists(path) else problems).append(f"{name}: {what} {'found' if os.path.exists(path) else 'MISSING'}: {path}")
        if any(s.startswith("mask_") or (s == "pr" and a3 != "era5") for s in steps):
            ok = check(f"{cof_root}cof_masks/{cof}_cofmasks_hp8_v1.zarr", f"{name}: COF store (Step 3 output, source of the masks and of tot_pr)")
            if ok:
                notes.append(f"{name}: COF store complete")
        if "pr" in steps and a3 == "era5":
            if check(defaults["imerg_6h_zarr"], f"{name}: IMERG 6-hourly store (precipitation)"):
                notes.append(f"{name}: IMERG 6-hourly store complete")
        if "link_env" in steps or "combine" in steps:
            src_dir = Path(f"{env_from}{a3}/single_vars") if "link_env" in steps else Path(info["single"])
            if not src_dir.is_dir():
                problems.append(f"{name}: environment stores folder missing: {src_dir}")
                continue
            excluded = set(info["env_exclude"])
            stores = sorted(n_ for n_ in os.listdir(src_dir) if n_.startswith("etc_2d_") and n_.endswith("_all_all.zarr")
                            and n_[len("etc_2d_"):-len("_all_all.zarr")] not in excluded)
            if len(stores) < defaults["min_env_stores"]:
                problems.append(f"{name}: only {len(stores)} environment stores in {src_dir}, at least {defaults['min_env_stores']} expected")
            good = sum(check(str(src_dir / s), f"{name}: environment store", n_pts) for s in stores)
            if good == len(stores) and stores:
                notes.append(f"{name}: {good} environment stores in {src_dir}: complete, {n_pts} storm points each")
    return problems, notes


# ---------------------------------------------------------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Run Analysis 3 (the ETC composites) for any set of sources, in dependency order",
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__.split("Examples")[1].split("Exit status")[0])
    p.add_argument("--data-root", required=False, help="root that receives the new products (etc_tracks/, etc_data/); use a test area first")
    p.add_argument("--cof-root", default=None, help="root of the Step 3 products that are read (cof_masks/; default: the registry's, the production tree)")
    p.add_argument("--env-from", default=None, help="folder with <source>/single_vars/ environment stores to link (default: the registry's)")
    p.add_argument("--sources", nargs="*", default=None, help="sources to run, names or aliases (default: all in the registry)")
    p.add_argument("--steps", nargs="*", default=None, help=f"only these steps ({', '.join(STEP_NAMES)}); the steps they need must exist")
    p.add_argument("--from", dest="from_step", default=None, help="only this step and the ones after it in the order " + " ".join(STEP_NAMES))
    p.add_argument("--registry", default=str(CONFIG_DIR / "config_etc_pipeline.yaml"))
    p.add_argument("--python", default=None, help="python that runs the steps (default: the registry's, the hackathon environment)")
    p.add_argument("--max-slots", type=int, default=None, help="CPU slots that may be busy at once (default: 80%% of the logical CPUs)")
    p.add_argument("--stagger-sec", type=float, default=10.0, help="minimum seconds between two step starts (default 10)")
    p.add_argument("--min-free-gb", type=float, default=60.0, help="do not start a step while less memory than this is available (default 60 GB)")
    p.add_argument("--resume", action="store_true", help="skip steps whose marker says success; re-run the rest")
    p.add_argument("--force", action="store_true", help="overwrite outputs that this runner did not create")
    p.add_argument("--dry-run", action="store_true", help="show the graph, the commands and the estimated schedule; run nothing")
    p.add_argument("--preflight-only", action="store_true", help="run the input checks and stop")
    p.add_argument("--skip-preflight", action="store_true", help="do not check the inputs first")
    p.add_argument("--step-args", nargs=2, action="append", metavar=("STEP", "ARGS"), default=[],
                   help='extra arguments appended to the command of one step (mask = all seven mask tasks), e.g. --step-args pr "--start_date 2020-03-01 '
                        '--end_date 2020-03-08" (repeatable); for tests, or a special run')
    p.add_argument("--log-dir", default=None, help="logs, status.json and resources.csv (default: <data-root>/pipeline_logs/<run id>)")
    return p.parse_args()


def main():
    args = parse_args()
    if not args.data_root:
        sys.exit("--data-root is required (a test area first; the production tree only when its outputs are meant to be replaced)")
    root = args.data_root.rstrip("/") + "/"
    logical = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    max_slots = args.max_slots or max(1, int(0.8 * logical))
    log_dir = args.log_dir or f"{root}pipeline_logs/{time.strftime('%Y%m%d_%H%M%S')}"

    defaults, sources = load_registry(args.registry)
    python = args.python or defaults["python"]
    cof_root = (args.cof_root or defaults["cof_root"]).rstrip("/") + "/"
    env_from = (args.env_from or defaults["env_from"]).rstrip("/") + "/"
    names = resolve_sources(args.sources, sources)
    steps = normalize_steps(args.steps) if args.steps else list(STEP_NAMES)
    if args.from_step:
        first = normalize_steps([args.from_step])[0]
        steps = [s for s in steps if STEP_NAMES.index(s) >= STEP_NAMES.index(first)]
    steps = [s for s in STEP_NAMES if s in steps]
    extra = {}
    for step_name, step_args in args.step_args:
        extra.setdefault(normalize_steps([step_name])[0], []).extend(shlex.split(step_args))

    all_tasks, infos = build_tasks(names, sources, defaults, root, python, cof_root, env_from, extra)
    tasks = select(all_tasks, steps)
    conflicts, missing, warnings = eng.prepare_plan(tasks, root, args.resume, args.force, lambda d: all_tasks[d].outputs if d in all_tasks else [])

    print(f"Data root: {root}\nCOF products read from: {cof_root}\nEnvironment stores linked from: {env_from}\nSources: {', '.join(names)}\n"
          f"Steps: {', '.join(steps)}\nPython: {python}\nSlot budget: {max_slots} of {logical} logical CPUs")
    for w in warnings:
        print(f"WARNING: {w}")
    if missing:
        print("\nA step needs an earlier step that is neither selected nor present:")
        for t, d in missing:
            print(f"  {t.label} needs {d[0]}/{d[1]}")
        print("Add the steps to the selection (--steps / --from) or run them first.")
        sys.exit(2)
    if conflicts:
        print("\nThese outputs exist and were not made by this runner (no marker); nothing is overwritten unless --force is given:")
        for t in conflicts:
            print(f"  {t.label}: {[o for o in t.outputs if os.path.exists(o)][0]}")
        sys.exit(2)

    if args.dry_run:
        times, span = eng.simulate(tasks, max_slots, args.stagger_sec / 60)
        print("\nGraph (start and end are estimates in minutes from the start; steps run in parallel where their needs allow):")
        for k in eng.topological_order(tasks):
            t = tasks[k]
            when = f"{times[k][0]:6.1f} -> {times[k][1]:6.1f}" if k in times else "   already done"
            print(f"  {t.label:44s} slots {t.slots:3d}  {when}   needs: {', '.join(f'{d[0]}/{d[1]}' for d in t.deps) or '-'}")
            print(f"      {' '.join(t.argv)}")
        print(f"\nEstimated wall time: {span:.0f} min ({span / 60:.2f} h) if the estimates hold and the node is not slowed down by the load.")
        return

    if not args.skip_preflight:
        problems, notes = preflight(tasks, infos, defaults, python, cof_root, env_from)
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

    git = eng.git_state()
    print(f"Code version: {git}\nLogs: {log_dir}  (status.json is updated every few seconds)\n", flush=True)
    good, text = eng.Runner(tasks, root, log_dir, max_slots, args.stagger_sec, args.min_free_gb, git).run()
    print("\n" + text)
    (Path(log_dir) / "summary.txt").write_text(text + "\n")
    sys.exit(0 if good else 1)


if __name__ == "__main__":
    main()
