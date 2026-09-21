"""Tests of the radii of scripts/calc_etc_spatial_stats.py.

The coordinates x and y of the ETC 2D stores are offsets in GRID POINTS (-80 ... 80 at 0.25 degrees). Until 2026-09-21 the script compared them with
the radius as it was given, so a radius of 10 "degrees" was a circle of 10 grid points (2.5 degrees). The radius is now in degrees:

  - the mask of radius 10 degrees at 0.25 degrees is a circle of 40 grid points (exact count, boundary points); 2.5 degrees is the old 10 points
  - a grid with different spacing along x and y (0.5 x 0.25) is converted with the right factor for each axis
  - lon_res and lat_res are read from the global attributes; a store without them, or with a non-positive value, is an error
  - the statistics use the degree circle: a spike 5 degrees from the centre is inside the 10-degree circle (it was outside before), one at 15 degrees is not;
    the fraction of a mask and the domain-mean precipitation are computed over the pixels of that circle
  - the output records radius_units = degrees and the resolution

Run:  python tests/test_etc_stats_radius.py      (or pytest tests/)
"""
import sys
from pathlib import Path

import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import calc_etc_spatial_stats as st  # noqa: E402


def grid(res=0.25, n=80):
    idx = np.arange(-n, n + 1)
    return xr.DataArray(idx, dims="x", coords={"x": idx}), xr.DataArray(idx, dims="y", coords={"y": idx})


def brute_count(radius_deg, lon_res, lat_res, n=80):
    return sum(1 for i in range(-n, n + 1) for j in range(-n, n + 1) if np.hypot(i * lon_res, j * lat_res) <= radius_deg)


def test_mask_is_a_circle_in_degrees():
    x, y = grid()
    m = st.create_circular_mask(x, y, 10.0, 0.25, 0.25)
    assert m.dims == ("y", "x") or set(m.dims) == {"x", "y"}
    assert int(m.sum()) == brute_count(10.0, 0.25, 0.25) and int(m.sum()) > 5000            # 40 grid points, about 5000 pixels
    at = lambda i, j: bool(m.sel(x=i, y=j))
    assert at(40, 0) and not at(41, 0) and at(0, 40) and not at(0, 41)
    assert at(28, 28) and not at(29, 29)                                                    # hypot 39.6 and 41.0
    # 2.5 degrees is what "10" meant before the fix: a circle of 10 grid points
    old = st.create_circular_mask(x, y, 2.5, 0.25, 0.25)
    assert (old == (np.hypot(x, y) <= 10)).all() and int(old.sum()) == brute_count(10.0, 1.0, 1.0)
    # a spike 5 degrees away (20 grid points) is inside the new circle, outside the old grid-point circle
    assert at(20, 0) and not bool((np.hypot(x, y) <= 10).sel(x=20, y=0))


def test_anisotropic_resolution():
    x, y = grid()
    m = st.create_circular_mask(x, y, 10.0, 0.5, 0.25)
    at = lambda i, j: bool(m.sel(x=i, y=j))
    assert at(20, 0) and not at(21, 0) and at(0, 40) and not at(0, 41)                     # 10 degrees along x is 20 points, along y 40 points
    assert at(14, 28) and not at(15, 28)                                                    # (7, 7) degrees inside, (7.5, 7) outside
    assert int(m.sum()) == brute_count(10.0, 0.5, 0.25)


def test_grid_resolution_from_attributes():
    ds = xr.Dataset(attrs={"lon_res": 0.5, "lat_res": 0.25})
    assert st.get_grid_resolution(ds) == (0.5, 0.25)
    for attrs, word in (({"lat_res": 0.25}, "lon_res"), ({"lon_res": 0.25}, "lat_res"), ({}, "lon_res")):
        try:
            st.get_grid_resolution(xr.Dataset(attrs=attrs)); raise AssertionError("missing attribute must fail")
        except ValueError as e:
            assert word in str(e)
    for attrs in ({"lon_res": 0.0, "lat_res": 0.25}, {"lon_res": 0.25, "lat_res": -0.25}):
        try:
            st.get_grid_resolution(xr.Dataset(attrs=attrs)); raise AssertionError("non-positive spacing must fail")
        except ValueError:
            pass


def make_ds():
    n = 80
    idx = np.arange(-n, n + 1)
    shape = (2, idx.size, idx.size)
    a = np.zeros(shape, dtype="float32"); a[0, n, n + 20] = 100.0; a[1, n, n + 60] = 100.0          # 5 degrees and 15 degrees from the centre, along x
    pr = np.ones(shape, dtype="float32")
    box = np.zeros(shape, dtype="float32"); box[:, n - 20:n + 21, n - 20:n + 21] = 5.0             # a 10 x 10 degree box (+-5 degrees), track id 5
    full = np.full(shape, 7.0, dtype="float32")
    ds = xr.Dataset(
        {"a": (("time", "y", "x"), a, {"units": "K"}), "pr": (("time", "y", "x"), pr, {"units": "mm h-1"}),
         "mcs_etc_overlap_mask": (("time", "y", "x"), box), "ar_etc_overlap_mask": (("time", "y", "x"), full),
         "overlap_flag": (("time",), np.array([1.0, 2.0]))},
        coords={"time": np.arange(2), "y": idx, "x": idx}, attrs={"lon_res": 0.25, "lat_res": 0.25, "radius": 20.0})
    return ds


def test_statistics_use_the_degree_circle():
    ds = make_ds()
    x, y = grid()
    circ10, circ15 = st.create_circular_mask(x, y, 10.0, 0.25, 0.25), st.create_circular_mask(x, y, 15.0, 0.25, 0.25)
    n10, n15 = int(circ10.sum()), int(circ15.sum())
    basic = st.compute_spatial_stats_basic(ds, ["a"], radius_deg=10)
    assert float(basic["a_max"][0]) == 100.0 and float(basic["a_max"][1]) == 0.0              # 5 degrees inside, 15 degrees outside
    assert abs(float(basic["a_mean"][0]) - 100.0 / n10) < 1e-6
    frac = st.compute_mask_fractional_area(ds, ["mcs_etc_overlap_mask", "ar_etc_overlap_mask"], radii=[10, 15])
    assert set(frac.data_vars) == {"mcs_etc_overlap_mask_frac_r10", "mcs_etc_overlap_mask_frac_r15", "ar_etc_overlap_mask_frac_r10", "ar_etc_overlap_mask_frac_r15"}
    box = ((np.abs(x) <= 20) & (np.abs(y) <= 20))
    assert abs(float(frac["mcs_etc_overlap_mask_frac_r10"][0]) - float((circ10 & box).sum()) / n10) < 1e-12
    assert abs(float(frac["mcs_etc_overlap_mask_frac_r15"][0]) - float((circ15 & box).sum()) / n15) < 1e-12
    assert float(frac["ar_etc_overlap_mask_frac_r10"][0]) == 1.0 and float(frac["ar_etc_overlap_mask_frac_r15"][0]) == 1.0
    prs = st.compute_feature_precipitation_stats(ds, radius_deg=10)
    assert abs(float(prs["pr_ar_mean"][0]) - 1.0) < 1e-6                                     # pr = 1 everywhere under a mask that covers everything
    assert abs(float(prs["pr_mcs_mean"][0]) - float((circ10 & box).sum()) / n10) < 1e-6      # domain-mean rain under the box mask


def test_output_records_the_units():
    out = st.compute_all_spatial_stats(make_ds(), basic_vars=["a"], mask_vars=["mcs_etc_overlap_mask"], basic_radius=10, mask_radii=[10, 15], pr_radius=10)
    assert out.attrs["radius_units"] == "degrees" and out.attrs["lon_res"] == 0.25 and out.attrs["lat_res"] == 0.25
    assert out.attrs["basic_radius_deg"] == 10 and list(out.attrs["mask_radii_deg"]) == [10, 15]
    no_res = make_ds(); no_res.attrs.pop("lon_res")
    try:
        st.compute_all_spatial_stats(no_res, basic_vars=["a"], mask_vars=[], basic_radius=10, mask_radii=[10], pr_radius=10); raise AssertionError("missing lon_res must fail")
    except ValueError:
        pass


if __name__ == "__main__":
    test_mask_is_a_circle_in_degrees(); test_anisotropic_resolution(); test_grid_resolution_from_attributes(); test_statistics_use_the_degree_circle(); test_output_records_the_units()
    print("test_etc_stats_radius: all checks passed")
