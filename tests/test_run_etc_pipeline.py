"""Tests of scripts/run_etc_pipeline.py (the graph, the commands, the protection) and of the engine change it needs in run_cof_pipeline.py.

  - 13 tasks per source: etc_cof, pr, seven mask_<variable>, link_env, combine, composites, stats; the dependencies of the graph
  - the commands: no stray tokens (an escaped newline of the Slurm-style text once became an argument of its own), outputs under the
    data root, the COF products read from the COF root, the expected share of missing storm points per source
  - the step selection ('mask' stands for the seven mask tasks), --step-args, source aliases
  - an output that exists without a marker is a conflict (protection) and no conflict with --force
  - --tracks-dir (every command reads the track file from that folder) and --extract-env (one env_<store> task per environment variable, named
    as the production stores, in place of link_env; the dependencies of combine; the step selection; link_env and env exclude each other)
  - --reuse-env: the same env_<store> tasks, each running subset_etc_env_store.py from the old store (--env-from) with the registry's track file as the
    old one and the --tracks-dir file as the new one; needs --tracks-dir; excludes --extract-env
  - the engine: the runner's data root is always COF_DATA_ROOT, whatever the caller's environment says, unless a task sets its own

Run:  python tests/test_run_etc_pipeline.py      (or pytest tests/)
"""
import os
import shutil
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import run_etc_pipeline as rp        # noqa: E402
import run_cof_pipeline as eng       # noqa: E402

ROOT, COF = "/tmp/etc_test_root/", "/tmp/cof_test_root/"


def graph(extra=None, names=None):
    rp.expected_points = lambda mc: 1000               # do not read the track files
    defaults, sources = rp.load_registry(str(REPO / "config" / "config_etc_pipeline.yaml"))
    names = names or list(sources)
    tasks, infos = rp.build_tasks(names, sources, defaults, ROOT, "/py/python", COF, "/env_from/", extra or {})
    return tasks, infos, defaults, sources


def opt(argv, name):
    return argv[argv.index(name) + 1]


def test_graph_and_commands():
    tasks, infos, defaults, sources = graph()
    steps = ["etc_cof", "pr", *[f"mask_{v}" for v in rp.MASK_VARS], "link_env", "combine", "composites", "stats"]
    assert len(sources) == 6 and len(tasks) == 6 * len(steps) == 78 and len(rp.MASK_VARS) == 7
    for name in sources:
        assert [k[1] for k in tasks if k[0] == name] == steps
        c = tasks[(name, "combine")]
        assert set(c.deps) == {(name, s) for s in ["etc_cof", "pr", *[f"mask_{v}" for v in rp.MASK_VARS], "link_env"]}
        assert tasks[(name, "composites")].deps == [(name, "combine")] and tasks[(name, "stats")].deps == [(name, "combine")]
        assert all(tasks[(name, s)].deps == [] for s in steps[:10])
    for t in tasks.values():
        assert all(a and a == a.strip() and "\n" not in a for a in t.argv), f"stray token in {t.label}: {t.argv}"
        assert all(o.startswith(ROOT) for o in t.outputs), t.outputs
    um = "um_glm_n2560_RAL3p3"
    pr, m = tasks[(um, "pr")], tasks[(um, "mask_mcs_ar_etc_overlap_mask")]
    assert opt(pr.argv, "--output_dir") == f"{ROOT}etc_data/{um}/single_vars/" and opt(pr.argv, "--variables") == "pr" and "--cof_mask" not in pr.argv
    assert opt(m.argv, "--variables") == "mcs_ar_etc_overlap_mask" and "--cof_mask" in m.argv
    assert pr.env["COF_DATA_ROOT"] == COF and m.env["COF_DATA_ROOT"] == COF and "COF_DATA_ROOT" not in tasks[(um, "combine")].env
    assert opt(pr.argv, "--max_missing_fraction") == "0.035" and opt(m.argv, "--max_missing_fraction") == "0.035"
    assert opt(tasks[("icon_d3hp003", "pr")].argv, "--max_missing_fraction") == "0.15"
    assert "--max_missing_fraction" not in tasks[("era5", "pr")].argv and "--structured_mesh" in tasks[("era5", "pr")].argv
    link = tasks[(um, "link_env")]
    assert opt(link.argv, "--expected-points") == "1000" and opt(link.argv, "--src-dir") == f"/env_from/{um}/single_vars"
    excl = link.argv[link.argv.index("--exclude") + 1:link.argv.index("--expected-points")]
    assert excl == ["pr", *rp.MASK_VARS]
    sc = tasks[("scream", "link_env")].argv
    assert sc[sc.index("--exclude") + 1:sc.index("--expected-points")] == ["pr", *rp.MASK_VARS, "ps"], "SCREAM's stale ps store is not linked"
    comb = tasks[(um, "combine")]
    assert opt(comb.argv, "--input_dir") == f"{ROOT}etc_data/{um}/single_vars/" and opt(comb.argv, "--etc-path") == f"{ROOT}etc_tracks/"
    assert comb.slots > tasks[(um, "pr")].slots and tasks[("era5", "pr")].est_min > tasks[(um, "pr")].est_min
    assert opt(tasks[(um, "etc_cof")].argv, "--cof_dir") == f"{COF}cof_masks/stats/" and opt(tasks[("era5", "etc_cof")].argv, "--source") == "era5"
    order = [k for k in eng.topological_order(tasks)]
    for k, t in tasks.items():
        assert all(order.index(d) < order.index(k) for d in t.deps)


def test_selection_extra_args_and_aliases():
    tasks, _, defaults, sources = graph(extra={"pr": ["--start_date", "2020-03-01"]})
    assert tasks[("scream", "pr")].argv[-2:] == ["--start_date", "2020-03-01"] and "--start_date" not in tasks[("scream", "combine")].argv
    assert {k[1] for k in rp.select(tasks, ["combine"])} == {"combine"}
    sel = rp.select(tasks, ["mask", "stats"])
    assert len([k for k in sel if k[0] == "scream"]) == 8
    assert rp.resolve_sources(["um", "UM", "obs", "icon"], sources) == ["um_glm_n2560_RAL3p3", "era5", "icon_d3hp003"]


ICON_ENV = {"tas", "huss", "ps", "psl", "uas", "vas", "prw",
            *[f"{v}_{l}hPa" for l in (850, 500) for v in ("ua", "va", "hus", "hur", "zg", "wa")]}       # the 19 environment stores of the March extraction
# SCREAM's production stores without pr, the seven masks and the stale ps (not extracted again)
SCREAM_ENV = {"hus_500hPa", "hus_850hPa", "huss", "omega500", "omega_500hPa", "omega850", "omega_850hPa", "psl", "rh500", "rh850", "tas", "ua500",
              "ua_500hPa", "ua850", "ua_850hPa", "uas", "uivt", "va500", "va_500hPa", "va850", "va_850hPa", "vas", "vivt", "zg500"}


def test_tracks_dir_and_extract_env():
    rp.expected_points = lambda mc: 1000
    defaults, sources = rp.load_registry(str(REPO / "config" / "config_etc_pipeline.yaml"))
    tasks, infos = rp.build_tasks(list(sources), sources, defaults, ROOT, "/py/python", COF, "/env_from/", {}, tracks_dir="/new/tracks/", extract_env=True)
    for name in sources:
        steps = [k[1] for k in tasks if k[0] == name]
        env = [s for s in steps if s.startswith("env_")]
        assert env
        assert "link_env" not in steps and steps[:2] == ["etc_cof", "pr"] and steps[-3:] == ["combine", "composites", "stats"]
        comb = tasks[(name, "combine")]
        assert set(comb.deps) == {(name, s) for s in ["etc_cof", "pr", *[f"mask_{v}" for v in rp.MASK_VARS], *env]}
        for s in env:                                                     # extracted from the catalog, one variable each, into the run's single_vars
            t = tasks[(name, s)]
            assert t.deps == [] and opt(t.argv, "--output_dir") == f"{ROOT}etc_data/{infos[name]['a3']}/single_vars/" and "--cof_mask" not in t.argv
            assert t.outputs == [f"{ROOT}etc_data/{infos[name]['a3']}/single_vars/etc_2d_{s[4:]}_all_all.zarr"]
            assert opt(t.argv, "--trackfile").startswith("/new/tracks/") and opt(t.argv, "--trackfile").endswith(".etc_stitched_nodes.filtered_out_tcs.txt")
        for s in ("pr", "mask_etc_ar_overlap_mask"):
            assert opt(tasks[(name, s)].argv, "--trackfile").startswith("/new/tracks/")
        assert opt(tasks[(name, "etc_cof")].argv, "--etc_dir") == "/new/tracks/"
        online = name in ("um_glm_n2560_RAL3p3", "casesm2_10km_nocumulus")        # online-only catalogs: at most 4 readers at a time (48 of 195 slots)
        assert all(t.slots == (48 if online else 4) for k, t in tasks.items() if k[0] == name and t.step.startswith("env_")) and tasks[(name, env[0])].est_min > 0
        if online:
            assert tasks[(name, env[0])].est_min == 50 and tasks[(name, env[0])].timeout_min == 240
    assert {t.step[4:] for k, t in tasks.items() if k[0] == "icon_d3hp003" and t.step.startswith("env_")} == ICON_ENV
    assert {t.step[4:] for k, t in tasks.items() if k[0] == "scream" and t.step.startswith("env_")} == SCREAM_ENV, "SCREAM's stale ps store is not extracted"
    assert len({t.step[4:] for k, t in tasks.items() if k[0] == "casesm2_10km_nocumulus" and t.step.startswith("env_")}) == 19
    # the default graph (no options) still reads the registry's track files and links
    tasks0, _ = rp.build_tasks(["icon_d3hp003"], sources, defaults, ROOT, "/py/python", COF, "/env_from/", {})
    assert ("icon_d3hp003", "link_env") in tasks0 and not any(k[1].startswith("env_") for k in tasks0)
    assert opt(tasks0[("icon_d3hp003", "pr")].argv, "--trackfile").startswith("/pscratch/sd/w/wcmca1/hackathon/etc_tracks/")
    # selection: 'env' stands for the env_<store> tasks and nothing else
    sel = rp.select(tasks, ["env"])
    assert sel and all(k[1].startswith("env_") for k in sel) and len([k for k in sel if k[0] == "icon_d3hp003"]) == 19
    assert {k[1] for k in rp.select(tasks, ["combine"])} == {"combine"}


def test_reuse_env_tasks():
    rp.expected_points = lambda mc: 1000
    defaults, sources = rp.load_registry(str(REPO / "config" / "config_etc_pipeline.yaml"))
    ext, _ = rp.build_tasks(["icon_d3hp003", "scream"], sources, defaults, ROOT, "/py/python", COF, "/env_from/", {}, tracks_dir="/new/tracks/", extract_env=True)
    reu, infos = rp.build_tasks(["icon_d3hp003", "scream"], sources, defaults, ROOT, "/py/python", COF, "/env_from/", {}, tracks_dir="/new/tracks/", reuse_env=True)
    assert set(reu) == set(ext), "the same steps as --extract-env, with the same outputs and dependencies"
    for k, t in reu.items():
        assert t.outputs == ext[k].outputs and t.deps == ext[k].deps and t.slots == ext[k].slots
    icon = [t for k, t in reu.items() if k[0] == "icon_d3hp003" and k[1].startswith("env_")]
    assert len(icon) == 19
    for t in icon:
        store = t.step[4:]
        assert t.argv[1].endswith("extract_environments/subset_etc_env_store.py") and "extract_etc_2d_vars.py" not in " ".join(t.argv)
        assert opt(t.argv, "--src-store") == f"/env_from/icon_d3hp003/single_vars/etc_2d_{store}_all_all.zarr"
        assert opt(t.argv, "--dst-store") == t.outputs[0] and t.outputs[0] == f"{ROOT}etc_data/icon_d3hp003/single_vars/etc_2d_{store}_all_all.zarr"
        assert opt(t.argv, "--new-track-file") == "/new/tracks/icon_d3hp003_hp8.etc_stitched_nodes.filtered_out_tcs.txt"
        assert opt(t.argv, "--old-track-file") == "/pscratch/sd/w/wcmca1/hackathon/etc_tracks/icon_d3hp003_hp8.etc_stitched_nodes.filtered_out_tcs.txt"
    assert {k[1][4:] for k in reu if k[0] == "scream" and k[1].startswith("env_")} == SCREAM_ENV, "SCREAM's stale ps store is not built"
    # the catalog of each variable's group (frame times, so that points without a frame become NaN as in the extraction)
    t3d = reu[("scream", "env_ua_850hPa")].argv
    assert opt(t3d, "--catalog-model") == "scream_ne120" and opt(t3d, "--catalog-url").endswith("main.yaml") and "--current-location" not in t3d
    assert opt(reu[("scream", "env_tas")].argv, "--catalog-model") == "scream_ne120_inst"
    assert opt(t3d, "--catalog-params") == '{"zoom": 8}'
    online, _ = rp.build_tasks(["casesm2_10km_nocumulus"], sources, defaults, ROOT, "/py/python", COF, "/env_from/", {}, tracks_dir="/new/tracks/", reuse_env=True)
    assert opt(online[("casesm2_10km_nocumulus", "env_wa_850hPa")].argv, "--current-location") == "online"
    assert infos["scream"]["old_points"] == 1000 and ("scream", "link_env") not in reu
    # pr, masks and etc_cof are unchanged: the new track file, the catalog / COF store
    assert opt(reu[("scream", "pr")].argv, "--trackfile").startswith("/new/tracks/") and "extract_etc_2d_vars.py" in " ".join(reu[("scream", "pr")].argv)


def test_reuse_env_option_rules():
    import subprocess
    base = [sys.executable, str(REPO / "scripts" / "run_etc_pipeline.py"), "--data-root", "/tmp/etc_test_root", "--dry-run"]
    bad = subprocess.run(base + ["--reuse-env"], capture_output=True, text=True)
    assert bad.returncode != 0 and "needs --tracks-dir" in (bad.stdout + bad.stderr)
    bad = subprocess.run(base + ["--reuse-env", "--extract-env", "--tracks-dir", "/tmp"], capture_output=True, text=True)
    assert bad.returncode != 0 and "exclude each other" in (bad.stdout + bad.stderr)
    bad = subprocess.run(base + ["--reuse-env", "--tracks-dir", "/tmp", "--steps", "link_env", "combine"], capture_output=True, text=True)
    assert bad.returncode != 0 and "exclude each other" in (bad.stdout + bad.stderr)


def test_link_env_and_env_exclude_each_other():
    import subprocess
    base = [sys.executable, str(REPO / "scripts" / "run_etc_pipeline.py"), "--data-root", "/tmp/etc_test_root", "--dry-run"]
    bad = subprocess.run(base + ["--extract-env", "--steps", "link_env", "combine"], capture_output=True, text=True)
    assert bad.returncode != 0 and "exclude each other" in (bad.stdout + bad.stderr)
    bad = subprocess.run(base + ["--steps", "env", "combine"], capture_output=True, text=True)
    assert bad.returncode != 0 and "exclude each other" in (bad.stdout + bad.stderr)
    bad = subprocess.run(base + ["--tracks-dir", "/no/such/folder"], capture_output=True, text=True)
    assert bad.returncode != 0 and "not a folder" in (bad.stdout + bad.stderr)


def test_protection_of_outputs_without_marker():
    tmp = tempfile.mkdtemp(prefix="etc_prot_")
    try:
        rp.expected_points = lambda mc: 1000
        defaults, sources = rp.load_registry(str(REPO / "config" / "config_etc_pipeline.yaml"))
        root = tmp + "/"
        tasks, _ = rp.build_tasks(["um_glm_n2560_RAL3p3"], sources, defaults, root, "/py/python", COF, "/env_from/", {})
        out = tasks[("um_glm_n2560_RAL3p3", "pr")].outputs[0]
        os.makedirs(out)                                              # exists, and no marker of the runner
        conflicts, missing, _ = eng.prepare_plan(tasks, root, False, False, lambda d: [])
        assert [t.label for t in conflicts] == ["um_glm_n2560_RAL3p3/pr"], [t.label for t in conflicts]
        conflicts, _, _ = eng.prepare_plan(tasks, root, False, True, lambda d: [])
        assert conflicts == [], "--force overwrites"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_engine_data_root_environment():
    """The runner's root is always COF_DATA_ROOT (even if the caller's environment has another), unless the task sets its own."""
    tmp = Path(tempfile.mkdtemp(prefix="etc_env_"))
    saved = os.environ.get("COF_DATA_ROOT")
    try:
        os.environ["COF_DATA_ROOT"] = "/caller/environment/"
        code = "import os,sys; open(sys.argv[1],'w').write(os.environ['COF_DATA_ROOT'])"
        tasks = {}
        for step, env in (("plain", {}), ("own", {"COF_DATA_ROOT": "/cof/products/"})):
            out = str(tmp / f"{step}.txt")
            tasks[("s", step)] = eng.Task("s", step, [sys.executable, "-c", code, out], env, 1, 0.1, 5.0, [out], [])
        root = str(tmp) + "/"
        good, _ = eng.Runner(tasks, root, tmp / "logs", 4, 0, 0, "test", poll_sec=0.05, sample_sec=1).run()
        assert good
        assert (tmp / "plain.txt").read_text() == root, "COF pipeline tasks must get the runner's root, not the caller's environment"
        assert (tmp / "own.txt").read_text() == "/cof/products/", "a task's own COF_DATA_ROOT wins"
    finally:
        if saved is None:
            os.environ.pop("COF_DATA_ROOT", None)
        else:
            os.environ["COF_DATA_ROOT"] = saved
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    test_graph_and_commands(); test_selection_extra_args_and_aliases(); test_tracks_dir_and_extract_env(); test_reuse_env_tasks(); test_reuse_env_option_rules()
    test_link_env_and_env_exclude_each_other()
    test_protection_of_outputs_without_marker(); test_engine_data_root_environment()
    print("test_run_etc_pipeline: all checks passed")
