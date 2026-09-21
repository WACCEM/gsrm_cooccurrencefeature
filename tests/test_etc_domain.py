"""Tests of the latitude domain of the ETC composites and statistics (src/etc_domain.py, create_etc_composites.py, calc_etc_spatial_stats.py).

The COF products (MCS masks, and with them the overlap flags and masks) exist only equatorward of 60 degrees, the ETC tracks are global. An ETC is used only
if at least min_coverage of its box (161 rows, +-20 degrees at 0.25 degrees) lies within |latitude| <= 60:

  - coverage of the box for centres at 0, 40, 48, 48.25, 60, 60.25 and 70 degrees, both hemispheres, NaN; the rounding to the grid
  - in_cof_domain at 0.8 (centre within 48 degrees), 0.5 (within 60 degrees, the "centroid" rule) and 1.0 (within 40); 0 keeps everything, also a NaN centre
  - the dataset version reads cof_lat (or storm_lat) and lat_res, and fails without them
  - the composites: the rule applies to every category (flag 3 and 'all'), in both hemispheres, with the counts and the mean of a marker field
    that identifies the points; the attributes record the rule; 0 switches it off
  - the statistics drop the points and keep the others in order

Run:  python tests/test_etc_domain.py      (or pytest tests/)
"""
import contextlib
import io
import sys
from pathlib import Path

import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from src import etc_domain as dom  # noqa: E402
import create_etc_composites as cc  # noqa: E402
import calc_etc_spatial_stats as st  # noqa: E402

ROWS = np.arange(-80, 81) * 0.25


def cov(lat):
    return float(dom.lat_box_coverage(np.array([lat]), ROWS, 0.25)[0])


def test_coverage_of_the_box():
    assert cov(0) == 1.0 and cov(40) == 1.0 and cov(-40) == 1.0                       # rows 20 ... 60, the 60.0 row is inside
    assert abs(cov(48) - 129 / 161) < 1e-12 and abs(cov(48.25) - 128 / 161) < 1e-12   # 80.1% and 79.5%
    assert abs(cov(60) - 81 / 161) < 1e-12 and abs(cov(60.25) - 80 / 161) < 1e-12     # 50.3% and 49.7%
    assert abs(cov(70) - 41 / 161) < 1e-12
    for lat in (30, 45, 48, 48.25, 55, 60, 60.25, 65):
        assert cov(lat) == cov(-lat), lat                                             # southern symmetry
    assert cov(np.nan) == 0.0
    assert cov(48.1) == cov(48.0) and cov(48.13) == cov(48.25)                        # the centre is rounded to the grid, as in the extraction


def test_selection_at_the_three_levels():
    lats = np.array([0.0, 40.0, 40.25, 48.0, 48.25, 55.0, 60.0, 60.25, 65.0, -48.0, -48.25, -60.0, -60.25, np.nan])
    sel = lambda c: dom.in_cof_domain(lats, ROWS, 0.25, 60.0, c)
    assert sel(0.8).tolist() == [True, True, True, True, False, False, False, False, False, True, False, False, False, False]
    assert sel(0.5).tolist() == [True, True, True, True, True, True, True, False, False, True, True, True, False, False]
    assert sel(1.0).tolist() == [True, True, False, False, False, False, False, False, False, False, False, False, False, False]
    assert sel(0.0).all() and sel(-1.0).all()                                         # switched off: every point, also a NaN centre
    # a different limit moves the boundary (the box of a centre at 40 degrees lies within 70 degrees up to row 60 only for limit >= 60)
    assert dom.in_cof_domain(np.array([50.0]), ROWS, 0.25, 70.0, 1.0)[0] and not dom.in_cof_domain(np.array([50.25]), ROWS, 0.25, 70.0, 1.0)[0]


def make_ds(lats, flags, name="cof_lat", attrs=None):
    n = len(lats)
    ny, nx = ROWS.size, 3
    pr = np.broadcast_to((np.arange(n) + 1.0)[:, None, None], (n, ny, nx)).astype("float32").copy()      # pr = index + 1 identifies the points
    ones = np.ones((n, ny, nx), dtype="float32")
    ds = xr.Dataset(
        {"pr": (("time", "y", "x"), pr, {"units": "mm h-1"}), "psl": (("time", "y", "x"), 100.0 * ones),
         "ar_mcs_etc_overlap_mask": (("time", "y", "x"), ones), "mcs_ar_etc_overlap_mask": (("time", "y", "x"), 0.0 * ones), "etc_mcs_ar_overlap_mask": (("time", "y", "x"), ones),
         "overlap_flag": (("time",), np.asarray(flags, dtype="float64")), name: (("time",), np.asarray(lats, dtype="float64")),
         "cof_lon": (("time",), np.zeros(n))},
        coords={"time": np.datetime64("2020-01-01") + np.arange(n) * np.timedelta64(6, "h"), "y": np.arange(-80, 81), "x": np.arange(-1, 2)},
        attrs={"lon_res": 0.25, "lat_res": 0.25, "radius": 20.0} if attrs is None else attrs)
    return ds


def test_dataset_version():
    ds = make_ds([30.0, 48.0, 48.25, -50.0], [3, 3, 3, 3])
    assert dom.store_in_cof_domain(ds, 60.0, 0.8).tolist() == [True, True, False, False]
    assert dom.store_in_cof_domain(ds, 60.0, 0.0).all()
    ds2 = make_ds([30.0, 48.0, 48.25, -50.0], [3, 3, 3, 3], name="storm_lat")            # falls back to storm_lat
    assert dom.store_in_cof_domain(ds2, 60.0, 0.8).tolist() == [True, True, False, False]
    for bad, word in ((make_ds([30.0], [3], name="other"), "cof_lat"), (make_ds([30.0], [3], attrs={"radius": 20.0}), "lon_res")):
        try:
            dom.store_in_cof_domain(bad, 60.0, 0.8); raise AssertionError("must fail")
        except ValueError as e:
            assert word in str(e)
    assert dom.store_in_cof_domain(make_ds([30.0], [3], attrs={}), 60.0, 0.0).all()      # rule off: no attributes needed


def run_composite(ds, flag, name, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return cc.create_composites(ds, overlap_category=cc.get_overlap_category(flag), overlap_flag=flag, overlap_name=name, **kw)


def test_composites_apply_the_rule_to_every_category():
    # index: 0 30N, 1 45N, 2 48N, 3 48.25N, 4 55N, 5 62N, 6 30S, 7 50S, 8 10N (tropics), 9 40N with flag 0;  pr = index + 1
    lats = [30, 45, 48, 48.25, 55, 62, -30, -50, 10, 40]
    flags = [3, 3, 3, 3, 3, 3, 3, 3, 3, 0]
    ds = make_ds(lats, flags)
    nh, sh, _ = run_composite(ds, 3, "3way", lat_limit=60.0, min_lat_coverage=0.8)
    assert nh.attrs["n_points"] == 3 and abs(float(nh["pr"].mean()) - 2.0) < 1e-6                          # indices 0, 1, 2
    assert sh.attrs["n_points"] == 1 and abs(float(sh["pr"].mean()) - 7.0) < 1e-6                          # index 6 (-50 has coverage 75%)
    assert nh.attrs["n_points_hemisphere_before_lat_rule"] == 7 and sh.attrs["n_points_hemisphere_before_lat_rule"] == 2
    assert nh.attrs["lat_limit"] == 60.0 and nh.attrs["min_lat_coverage"] == 0.8
    nh, sh, _ = run_composite(ds, 3, "3way", lat_limit=60.0, min_lat_coverage=0.5)                          # centroid rule
    assert nh.attrs["n_points"] == 5 and abs(float(nh["pr"].mean()) - 3.0) < 1e-6
    assert sh.attrs["n_points"] == 2 and abs(float(sh["pr"].mean()) - 7.5) < 1e-6
    nh, sh, _ = run_composite(ds, 3, "3way", lat_limit=60.0, min_lat_coverage=0.0)                          # off
    assert nh.attrs["n_points"] == 6 and abs(float(nh["pr"].mean()) - 3.5) < 1e-6 and sh.attrs["n_points"] == 2
    nh, sh, _ = run_composite(ds, -1, "all", lat_limit=60.0, min_lat_coverage=0.8)                          # 'all' is filtered too
    assert nh.attrs["n_points"] == 4 and abs(float(nh["pr"].mean()) - 4.0) < 1e-6                          # indices 0, 1, 2 and 9
    nh, sh, _ = run_composite(ds, 0, "isolated", lat_limit=60.0, min_lat_coverage=0.8)
    assert nh.attrs["n_points"] == 1 and abs(float(nh["pr"].mean()) - 10.0) < 1e-6


def test_statistics_drop_the_points():
    lats = [30.0, 58.0, 60.0, 60.25, 65.0, -59.0, -61.0]
    ds = make_ds(lats, [0, 1, 2, 3, 0, 1, 2])
    kept, n_before, n_kept = st.apply_lat_domain(ds, 60.0, 0.5)                          # centroid rule
    assert (n_before, n_kept) == (7, 4) and kept.cof_lat.values.tolist() == [30.0, 58.0, 60.0, -59.0]
    assert kept.overlap_flag.values.tolist() == [0, 1, 2, 1] and kept.sizes["time"] == 4
    kept, n_before, n_kept = st.apply_lat_domain(ds, 60.0, 0.8)                          # the composites' sample
    assert n_kept == 1 and kept.cof_lat.values.tolist() == [30.0]
    kept, n_before, n_kept = st.apply_lat_domain(ds, 60.0, 0.0)                          # off: the same object comes back
    assert n_kept == 7 and kept is ds


if __name__ == "__main__":
    test_coverage_of_the_box(); test_selection_at_the_three_levels(); test_dataset_version(); test_composites_apply_the_rule_to_every_category(); test_statistics_drop_the_points()
    print("test_etc_domain: all checks passed")
