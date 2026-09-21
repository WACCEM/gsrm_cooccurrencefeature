"""Tests of scripts/run_etc_pipeline.py (the graph, the commands, the protection) and of the engine change it needs in run_cof_pipeline.py.

  - 13 tasks per source: etc_cof, pr, seven mask_<variable>, link_env, combine, composites, stats; the dependencies of the graph
  - the commands: no stray tokens (an escaped newline of the Slurm-style text once became an argument of its own), outputs under the
    data root, the COF products read from the COF root, the expected share of missing storm points per source
  - the step selection ('mask' stands for the seven mask tasks), --step-args, source aliases
  - an output that exists without a marker is a conflict (protection) and no conflict with --force
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
    test_graph_and_commands(); test_selection_extra_args_and_aliases(); test_protection_of_outputs_without_marker(); test_engine_data_root_environment()
    print("test_run_etc_pipeline: all checks passed")
