"""Tests of the thresholds version in scripts/run_cof_pipeline.py and config/config_pipeline.yaml (2026-09-23):

  - IMERG's registry entry makes the thresholds from the 11-year (2014-2024) 6-hourly store with --version v1_2014_2024, the
    thresholds output has that version in its name, and the attribution is given that same file through --threshold_file
  - a source without a thresholds version (the five models) is unchanged: the default file name, no --version, no --threshold_file
  - a version in the registry works for any source, and the file the attribution reads is always the file the thresholds step writes
  - arguments given with --step-args come after the registry's, so they override them (the last occurrence wins in argparse)
  - the attribution's own default threshold file (threshold_version in config_sources.yaml) is the file the runner's thresholds
    step writes, for every source: the two settings cannot drift apart

Run:  python tests/test_run_cof_pipeline_thresholds.py      (or pytest tests/)
"""
import copy
import os
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import run_cof_pipeline as rp  # noqa: E402
import calc_stormtype_extreme_precip_spatial as sp  # noqa: E402

ROOT = "/data/root/"
IMERG_11YR_STORE = "/pscratch/sd/w/wcmca1/GPM/healpix/IMERG_V7_6H_zoom8_20140101_20241231.zarr"


def registry():
    return rp.load_registry(REPO / "config" / "config_pipeline.yaml")


def build(names, steps, sources, defaults, extra_args=None):
    return rp.build_tasks(names, steps, sources, defaults, ROOT, "python", extra_args=extra_args)


def last(argv, flag):
    """Value after the last occurrence of an option: what argparse uses when it is given more than once."""
    return argv[len(argv) - 1 - argv[::-1].index(flag) + 1]


def test_imerg_thresholds_come_from_the_11_year_store_and_the_attribution_reads_that_file():
    defaults, sources = registry()
    tasks = build(["IMERGv7"], ["thresholds", "attribution"], sources, defaults)
    thr, att = tasks[("IMERGv7", "thresholds")], tasks[("IMERGv7", "attribution")]

    assert last(thr.argv, "--input_zarr") == IMERG_11YR_STORE and last(thr.argv, "--input_var") == "precipitation"
    assert last(thr.argv, "--version") == "v1_2014_2024"
    assert thr.outputs == [f"{ROOT}extreme_precip/IMERGv7_precip_percentiles_6h_hp8_v1_2014_2024.nc"]
    assert thr.deps == []                                      # IMERG's thresholds read their own store, no earlier step
    # the attribution reads the file that the thresholds step writes, and needs that step
    assert last(att.argv, "--threshold_file") == thr.outputs[0]
    assert ("IMERGv7", "thresholds") in att.deps and ("IMERGv7", "s3") in att.deps
    assert att.outputs == [f"{ROOT}extreme_precip/IMERGv7_stormtype_spatial_p90.nc", f"{ROOT}extreme_precip/IMERGv7_stormtype_spatial_p95.nc"]


def test_the_models_are_unchanged():
    defaults, sources = registry()
    models = [n for n in sources if n != "IMERGv7"]
    assert len(models) == 5
    tasks = build(models, ["thresholds", "attribution"], sources, defaults)
    for name in models:
        thr, att = tasks[(name, "thresholds")], tasks[(name, "attribution")]
        assert "--version" not in thr.argv and "--threshold_file" not in att.argv, name
        assert thr.outputs == [f"{ROOT}extreme_precip/{sources[name]['source_name']}_precip_percentiles_6h_hp8_v1.nc"], name
        assert (name, "s1") in thr.deps, name                   # they read Step 1's tot_pr


def test_a_version_in_the_registry_works_for_any_source():
    defaults, sources = registry()
    sources = copy.deepcopy(sources)
    sources["scream"]["thresholds"]["version"] = "v9"
    tasks = build(["scream"], ["thresholds", "attribution"], sources, defaults)
    thr, att = tasks[("scream", "thresholds")], tasks[("scream", "attribution")]
    assert last(thr.argv, "--version") == "v9"
    assert thr.outputs == [f"{ROOT}extreme_precip/scream_precip_percentiles_6h_hp8_v9.nc"]
    assert last(att.argv, "--threshold_file") == thr.outputs[0]
    assert rp.thresholds_file(sources["scream"], ROOT) == thr.outputs[0]


def test_the_attribution_default_threshold_file_is_the_file_the_runner_makes():
    defaults, sources = registry()
    cfg_sources = yaml.safe_load(open(REPO / "config" / "config_sources.yaml"))
    for name, src in sources.items():
        made_by_runner = os.path.normpath(rp.thresholds_file(src, ROOT))
        read_by_default = os.path.normpath(sp.default_threshold_file(cfg_sources[src["config_key"]], ROOT))
        assert read_by_default == made_by_runner, name
    # and IMERG's is the 2014-2024 file
    assert cfg_sources["IR_IMERG"]["threshold_version"] == "v1_2014_2024"
    assert all("threshold_version" not in cfg_sources[s["config_key"]] for n, s in sources.items() if n != "IMERGv7")


def test_step_args_override_the_registry():
    defaults, sources = registry()
    extra = {"thresholds": ["--input_zarr", "/other/store.zarr", "--version", "v20yr"],
             "attribution": ["--threshold_file", "/other/thresholds.nc"]}
    tasks = build(["IMERGv7"], ["thresholds", "attribution"], sources, defaults, extra_args=extra)
    thr, att = tasks[("IMERGv7", "thresholds")], tasks[("IMERGv7", "attribution")]
    assert last(thr.argv, "--input_zarr") == "/other/store.zarr" and last(thr.argv, "--version") == "v20yr"
    assert last(att.argv, "--threshold_file") == "/other/thresholds.nc"


if __name__ == "__main__":
    test_imerg_thresholds_come_from_the_11_year_store_and_the_attribution_reads_that_file()
    test_the_models_are_unchanged()
    test_a_version_in_the_registry_works_for_any_source()
    test_the_attribution_default_threshold_file_is_the_file_the_runner_makes()
    test_step_args_override_the_registry()
    print("test_run_cof_pipeline_thresholds: all checks passed")
