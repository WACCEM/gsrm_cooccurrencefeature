"""Tests of scripts/compare_processing_rounds.py (exact comparison of the products of two processing rounds).

  - identical roots: nothing differs, exit status 0 with --expect identical
  - a store with changed cells: the variable, the number of differing cell-times and the two footprints are reported, an unchanged variable is not
  - a netCDF variable and a parquet column that changed, a product that exists in one root only, exit status 1 with --expect identical

Run:  python tests/test_compare_processing_rounds.py      (or pytest tests/)
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

REPO = Path(__file__).resolve().parent.parent
SCRIPT = str(REPO / "scripts" / "compare_processing_rounds.py")
sys.path.insert(0, str(REPO / "scripts"))
import compare_processing_rounds as cmp   # noqa: E402

SRC = "srcx"


def make_root(root, tc_extra=0, pr_shift=0.0, etc_col=0):
    (root / "cof_masks" / "stats" / "monthly").mkdir(parents=True)
    (root / "extreme_precip").mkdir()
    n_t, n_c = 6, 20
    tc = np.zeros((n_t, n_c), dtype="float32"); tc[1, :4] = 7; tc[2, 5:6 + tc_extra] = 3            # tc_extra more cells in the new root
    ar = np.zeros((n_t, n_c), dtype="float32"); ar[0, :3] = 1
    ds = xr.Dataset({"tc_mask": (("time", "cell"), tc), "ar_mask": (("time", "cell"), ar)},
                    coords={"time": np.arange(n_t).astype("datetime64[h]").astype("datetime64[ns]"), "cell": np.arange(n_c)})
    ds.to_zarr(root / "cof_masks" / f"{SRC}_cofmasks_hp8_v1.zarr", mode="w")
    xr.Dataset({"pr_p95": (("cell",), np.linspace(0, 1, n_c) + pr_shift), "pr_p90": (("cell",), np.linspace(0, 2, n_c))}).to_netcdf(
        root / "extreme_precip" / f"{SRC}_precip_percentiles_6h_hp8_v1.nc")
    pd.DataFrame({"etc_track": [1, 2, 3], "overlap_flag": [0, 1, 0] if etc_col == 0 else [0, 1, 1]}).to_parquet(root / "cof_masks" / "stats" / f"{SRC}_etc_overlap_tracking.parquet")


def test_identical_and_different_roots():
    tmp = Path(tempfile.mkdtemp(prefix="cmpr_"))
    try:
        a, b, c = tmp / "a", tmp / "b", tmp / "c"
        make_root(a); make_root(b); make_root(c, tc_extra=2, pr_shift=0.5, etc_col=1)
        lines = []
        n = cmp.run(str(a) + "/", str(b) + "/", [SRC], workers=2, out=lines.append)
        assert n == 0 and all("IDENTICAL" in l for l in lines if l.strip().startswith(("IDENTICAL", "DIFFERS"))) and len([l for l in lines if "IDENTICAL" in l]) == 3, lines
        lines = []
        n = cmp.run(str(a) + "/", str(c) + "/", [SRC], workers=2, out=lines.append)
        text = "\n".join(lines)
        assert n == 3, text                                                          # the store, the netCDF file and the parquet file differ
        assert "tc_mask" in text and "cell-times differing" in text
        store_line = [l for l in lines if "tc_mask" in l][0]
        assert "2" in store_line.split("cell-times differing")[1].split("footprint")[0] and "footprint old            5 new            7" in store_line, store_line
        assert not any("ar_mask" in l for l in lines), "an unchanged variable is not listed"
        assert "pr_p95" in text and "overlap_flag" in text
        # a product that exists in one root only
        (a / "cof_masks" / "stats" / f"{SRC}_mcs_cof_flags.parquet").write_bytes((a / "cof_masks" / "stats" / f"{SRC}_etc_overlap_tracking.parquet").read_bytes())
        lines = []
        n = cmp.run(str(a) + "/", str(b) + "/", [SRC], workers=2, out=lines.append)
        assert n == 1 and any("MISSING" in l and "only in the old root" in l for l in lines), lines
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_exit_status_with_expect_identical():
    tmp = Path(tempfile.mkdtemp(prefix="cmpr_"))
    try:
        a, b, c = tmp / "a", tmp / "b", tmp / "c"
        make_root(a); make_root(b); make_root(c, tc_extra=1)
        base = [sys.executable, SCRIPT, "--old-root", str(a), "--sources", SRC, "--workers", "2"]
        r = subprocess.run(base + ["--new-root", str(b), "--expect", "identical"], capture_output=True, text=True)
        assert r.returncode == 0 and "RESULT: 0 product(s) differ" in r.stdout, r.stdout + r.stderr
        r = subprocess.run(base + ["--new-root", str(c), "--expect", "identical"], capture_output=True, text=True)
        assert r.returncode == 1 and "DIFFERS" in r.stdout
        r = subprocess.run(base + ["--new-root", str(c)], capture_output=True, text=True)
        assert r.returncode == 0, "without --expect the differences are only reported"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    test_identical_and_different_roots(); test_exit_status_with_expect_identical()
    print("test_compare_processing_rounds: all checks passed")
