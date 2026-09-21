"""Tests of extract_environments/combine_etc_2d_vars.py: unit-aware scaling, alignment of the single-variable stores, strict loading.

  - a native precipitation flux is scaled to mm/h with the source's factor (as before)
  - a `pr` that is already in mm h-1 (taken from Step 1's tot_pr) is NOT scaled again
  - a `pr` that is off by orders of magnitude after the standardization is an error (magnitude check)
  - stores with different storm points (count, times, IDs) are an error, not a warning followed by a merge by position
  - a variable that cannot be loaded (hollow store) is an error unless --allow-missing-variables

Run:  python tests/test_etc_combine.py      (or pytest tests/)
"""
import contextlib
import io
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "extract_environments"))
sys.path.insert(0, str(REPO / "scripts"))
import combine_etc_2d_vars as cb  # noqa: E402

N, NY, NX = 6, 5, 5
TIMES = np.datetime64("2020-01-01T00:00", "ns") + np.arange(N) * np.timedelta64(6, "h")


def make_store(path, var, value, units, times=TIMES, ids=None, n=N):
    ids = np.arange(1, n + 1) if ids is None else ids
    ds = xr.Dataset(
        {var: (("time", "y", "x"), np.full((n, NY, NX), value, dtype="float32"), {"units": units}),
         "storm_id": (("time",), ids), "storm_lat": (("time",), np.linspace(30, 40, n)), "storm_lon": (("time",), np.linspace(10, 20, n)),
         "grid_id": (("time",), np.arange(n))},
        coords={"time": times[:n], "y": np.arange(-2, 3), "x": np.arange(-2, 3)})
    ds.to_zarr(str(path), mode="w")
    return str(path)


def combine(files, out, source="scream", **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return cb.combine_zarr_files(files, str(out), "all_all", source=source, add_cof_data=False, **kw)


def mean_pr(out):
    return float(xr.open_zarr(str(out))["pr"].mean())


def test_all():
    tmp = Path(tempfile.mkdtemp(prefix="cmb_"))
    try:
        psl = make_store(tmp / "etc_2d_psl_all_all.zarr", "psl", 101000.0, "Pa")
        # 1) native SCREAM flux (m s-1): scaled by 3.6e6 to mm/h, as before
        pr_native = make_store(tmp / "native.zarr", "pr", 2.78e-8, "m s^-1")
        out = tmp / "c1.zarr"; combine([("pr", pr_native), ("psl", psl)], out)
        assert abs(mean_pr(out) - 0.1) < 0.002, mean_pr(out)
        assert xr.open_zarr(str(out))["pr"].attrs["scale_factor_applied"] == 3600000.0
        # 2) NICAM native flux (kg m-2 s-1): x3600
        pr_kg = make_store(tmp / "kg.zarr", "pr", 2.78e-5, "kg m-2 s-1")
        out = tmp / "c2.zarr"; combine([("pr", pr_kg), ("psl", psl)], out, source="nicam_gl11")
        assert abs(mean_pr(out) - 0.1) < 0.002
        # 3) pr from tot_pr (already mm h-1): NOT scaled again
        pr_tot = make_store(tmp / "tot.zarr", "pr", 0.1, "mm h-1")
        out = tmp / "c3.zarr"; combine([("pr", pr_tot), ("psl", psl)], out)
        a = xr.open_zarr(str(out))["pr"].attrs
        assert abs(mean_pr(out) - 0.1) < 1e-6, f"a mm h-1 pr must not be scaled again (mean {mean_pr(out)})"
        assert a["scale_factor_applied"] == 1.0 and a["units"] == "mm h-1"
        # 4) implausible magnitude after the standardization is an error
        pr_bad = make_store(tmp / "bad.zarr", "pr", 360.0, "mm h-1")
        try:
            combine([("pr", pr_bad), ("psl", psl)], tmp / "c4.zarr"); raise AssertionError("an implausible pr must fail")
        except ValueError as e:
            assert "outside" in str(e)
        # 5) storm points that differ: count, time, ID
        for label, kw in (("count", dict(n=N - 1)), ("times", dict(times=TIMES + np.timedelta64(6, "h"))), ("ids", dict(ids=np.arange(11, 11 + N)))):
            other = make_store(tmp / f"mis_{label}.zarr", "pr", 0.1, "mm h-1", **kw)
            try:
                combine([("pr", other), ("psl", psl)], tmp / f"c5_{label}.zarr"); raise AssertionError(f"different {label} must fail")
            except cb.AlignmentError:
                pass
        # 6) a hollow store: error by default, skipped only on request
        hollow = tmp / "hollow.zarr"; (hollow / "pr").mkdir(parents=True)
        try:
            combine([("pr", pr_tot), ("psl", psl), ("hus", str(hollow))], tmp / "c6.zarr"); raise AssertionError("a hollow store must fail")
        except RuntimeError as e:
            assert "hus" in str(e)
        out = tmp / "c6b.zarr"; combine([("pr", pr_tot), ("psl", psl), ("hus", str(hollow))], out, allow_missing_variables=True)
        assert set(xr.open_zarr(str(out)).data_vars) >= {"pr", "psl"} and "hus" not in xr.open_zarr(str(out)).data_vars
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    test_all()
    print("test_etc_combine: all checks passed")
