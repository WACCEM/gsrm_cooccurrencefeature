"""Tests of the per-calendar-year extreme precipitation amounts by storm type added to
scripts/calc_stormtype_extreme_precip_spatial.py (2026-09-23), after the per-year percentiles of calc_extreme_precip_thresholds.py:

  - qualifying_years: the coverage gate counts distinct calendar days (not time steps), keeps a year with exactly the
    minimum number of days and drops one with a day fewer, works for cftime calendars, and returns nothing below 2 years
  - process_timeseries_dask: the whole-record variables are identical (values, dtypes, attributes) with and without the
    per-year option, and equal a plain step-by-step float32 accumulation in time order, so the per-year accumulators sit next
    to the whole-record ones without touching them
  - the per-year variables: exactly the year coordinate plus total_extreme_count_annual, total_extreme_precip_annual and
    {type}_precip_annual for the 13 types; each year equals a direct sum over that year's time steps; the storm types add up to
    the year's total; the years add up to the whole record
  - save_spatial_results: year is int32, the amounts are compressed float32, the counts int32, the existing variables keep
    their encoding, the extra global attributes come after the standard ones
  - the command line: defaults (360 days, per-year output on, default threshold file), --threshold_file, --no_annual and
    --min_year_coverage_days, end to end through main() on a tiny synthetic data root
  - the default threshold file follows the source's threshold_version in config_sources.yaml (v1 when it has none)

Run:  python tests/test_stormtype_annual.py      (or pytest tests/)
"""
import os
import sys
import tempfile
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import calc_stormtype_extreme_precip_spatial as sp  # noqa: E402

TYPE_KEYS = ['mcs_isolated', 'ar_isolated', 'etc_isolated', 'tc', 'mcs_ar_2way', 'mcs_etc_2way', 'ar_etc_2way',
             'mcs_ar_etc_3way', 'dc', 'nd', 'st', 'dz', 'unassigned']
MASK_VARS = ['mcs_ar_overlap_mask', 'ar_mcs_overlap_mask', 'mcs_etc_overlap_mask', 'etc_mcs_overlap_mask',
             'ar_etc_overlap_mask', 'etc_ar_overlap_mask', 'mcs_ar_etc_overlap_mask', 'ar_mcs_etc_overlap_mask',
             'etc_mcs_ar_overlap_mask', 'mcs_isolated_mask', 'ar_isolated_mask', 'etc_isolated_mask', 'tc_mask']
ANNUAL_VARS = {'total_extreme_count_annual', 'total_extreme_precip_annual', *[f'{k}_precip_annual' for k in TYPE_KEYS]}
YEARS = [2019, 2020, 2021]


def monthly_times(years=YEARS):
    """One time step at the start of every month: enough to put 12 steps in each calendar year cheaply."""
    return pd.date_range(f"{years[0]}-01-01", f"{years[-1]}-12-01", freq="MS")


def make_case(times, n_cells=40, seed=0, quantile=0.75):
    """Synthetic (time, cell) precipitation, the storm masks that process_single_timestep reads (about 8% of the cells per
    mask, so all the categories, the overlaps and the unassigned cells occur) and per-cell thresholds."""
    rng = np.random.default_rng(seed)
    n = len(times)
    coords = {"time": times, "cell": np.arange(n_cells)}
    pr = xr.DataArray(rng.gamma(2.0, 0.5, size=(n, n_cells)).astype("float32"), dims=("time", "cell"),
                      coords=coords, name="tot_pr")
    masks = {name: (("time", "cell"), (rng.random((n, n_cells)) < 0.08).astype("float32")) for name in MASK_VARS}
    masks["cloud_types"] = (("time", "cell"), rng.integers(0, 5, size=(n, n_cells)).astype("float32"))
    ds = xr.Dataset(masks, coords=coords)
    threshold = np.quantile(pr.values, quantile, axis=0).astype("float32")
    return pr, ds, threshold


def run(pr, thr, ds, **kw):
    kw.setdefault("percentile_name", "P95")
    kw.setdefault("n_workers", 2)
    kw.setdefault("batch_size", 7)  # not a divisor of the number of steps, so a batch boundary falls inside a year
    return sp.process_timeseries_dask(pr, thr, ds, **kw)


def time_axis(times):
    return xr.DataArray(np.zeros(len(times)), dims="time", coords={"time": times})["time"]


# ---------------------------------------------------------------------------------------------------------------------
def test_qualifying_years_counts_distinct_days_and_needs_two_years():
    full = pd.date_range("2019-01-01", "2021-12-31 18:00", freq="6h")
    assert sp.qualifying_years(time_axis(full), 360) == [2019, 2020, 2021]
    # counts days, not time steps: 3-hourly and 6-hourly data give the same answer
    assert sp.qualifying_years(time_axis(pd.date_range("2019-01-01", "2021-12-31 21:00", freq="3h")), 360) == [2019, 2020, 2021]

    # a year with 60 days of data is dropped, not just given fewer samples
    short = full[~((full.year == 2021) & (full.dayofyear > 60))]
    assert sp.qualifying_years(time_axis(short), 360) == [2019, 2020]

    # exactly the minimum number of distinct days qualifies (>=), one day fewer does not; here 2019 has 365 - 5 = 360 days
    gap = full[~((full.year == 2019) & (full.dayofyear > 360))]
    assert sp.qualifying_years(time_axis(gap), 360) == [2019, 2020, 2021]
    gap = full[~((full.year == 2019) & (full.dayofyear > 359))]
    assert sp.qualifying_years(time_axis(gap), 360) == [2020, 2021]

    # below 2 qualifying years nothing is returned, even for a complete year
    assert sp.qualifying_years(time_axis(pd.date_range("2020-01-01", "2020-12-31 18:00", freq="6h")), 360) == []
    assert sp.qualifying_years(time_axis(short), 360, min_years=3) == []

    # a record that crosses a calendar boundary (about one year of a model) has no year with enough days
    crossing = pd.date_range("2019-08-01", "2020-08-31 18:00", freq="6h")
    assert sp.qualifying_years(time_axis(crossing), 360) == []


def test_qualifying_years_with_a_cftime_calendar():
    import cftime  # noqa: F401
    t = xr.cftime_range("2019-01-01", "2021-12-30 18:00", freq="6h", calendar="noleap")
    assert sp.qualifying_years(time_axis(t), 360) == [2019, 2020, 2021]


# ---------------------------------------------------------------------------------------------------------------------
def test_whole_record_variables_are_identical_with_and_without_the_annual_option():
    pr, ds, thr = make_case(monthly_times())
    base = run(pr, thr, ds)
    with_annual = run(pr, thr, ds, annual_years=YEARS)

    assert set(with_annual.variables) - set(base.variables) == ANNUAL_VARS | {"year"}
    assert not (set(base.variables) - set(with_annual.variables))
    # values, dtypes, dimensions, coordinates and attributes of every existing variable, and the global attributes
    xr.testing.assert_identical(with_annual[list(base.data_vars)], base)
    for name in base.data_vars:
        assert with_annual[name].dtype == base[name].dtype, name
    # no per-year request (None or empty) adds nothing
    for none in (None, []):
        xr.testing.assert_identical(run(pr, thr, ds, annual_years=none), base)


def test_whole_record_values_equal_a_plain_step_by_step_accumulation():
    """The float32 accumulators, in time order, as before the per-year addition: guards their order and dtype."""
    pr, ds, thr = make_case(monthly_times())
    out = run(pr, thr, ds, annual_years=YEARS)

    template = pr.isel(time=0).compute()
    counts = {k: xr.zeros_like(template) for k in ["total_extreme"] + TYPE_KEYS}
    precip = {k: xr.zeros_like(template) for k in ["total_extreme"] + TYPE_KEYS}
    for t in range(pr.sizes["time"]):
        r = sp.process_single_timestep(pr.isel(time=t), ds.isel(time=t), thr, True)
        counts["total_extreme"] += r["extreme_mask"]
        precip["total_extreme"] += r["total_extreme_pr"]
        for k in TYPE_KEYS:
            counts[k] += r[k]
            precip[k] += r[f"{k}_pr"]

    assert out["total_extreme_precip"].dtype == np.float32
    assert np.array_equal(out["total_extreme_count"].values, counts["total_extreme"].values)
    assert np.array_equal(out["total_extreme_precip"].values, precip["total_extreme"].values)
    for k in TYPE_KEYS:
        assert np.array_equal(out[f"{k}_count"].values, counts[k].values), k
        frac = xr.where(precip["total_extreme"] > 0, precip[k] / precip["total_extreme"], 0.0)
        assert np.array_equal(out[f"{k}_frac"].values, frac.values), k


def test_annual_totals_match_direct_sums_and_close_the_whole_record():
    times = monthly_times()
    pr, ds, thr = make_case(times)
    out = run(pr, thr, ds, annual_years=YEARS)

    assert list(out["year"].values) == YEARS and out["year"].dtype == np.int32
    for name in ANNUAL_VARS:
        assert out[name].dims == ("year", "cell") and out[name].dtype == np.float32, name

    step_year = np.asarray(times.year)
    for i, y in enumerate(YEARS):
        steps = np.flatnonzero(step_year == y)
        extreme = pr.values[steps] > thr
        assert np.array_equal(out["total_extreme_count_annual"].isel(year=i).values, extreme.sum(0))
        direct = np.where(extreme, pr.values[steps].astype("float64"), 0.0).sum(0)
        assert np.allclose(out["total_extreme_precip_annual"].isel(year=i).values, direct, rtol=1e-5, atol=1e-6)
        # per storm type: the type assignment of each step from process_single_timestep, summed over the year's steps
        per_step = [sp.process_single_timestep(pr.isel(time=t), ds.isel(time=t), thr, True) for t in steps]
        for k in TYPE_KEYS:
            want = sum(r[f"{k}_pr"].values.astype("float64") for r in per_step)
            assert np.allclose(out[f"{k}_precip_annual"].isel(year=i).values, want, rtol=1e-5, atol=1e-6), (y, k)

    # every extreme cell is assigned to exactly one type (unassigned included), so the types add up to the year's total
    by_type = sum(out[f"{k}_precip_annual"] for k in TYPE_KEYS)
    assert np.allclose(by_type.values, out["total_extreme_precip_annual"].values, rtol=1e-5, atol=1e-6)
    # ... and the years add up to the whole record (counts exactly, amounts up to the float32 order of summation)
    assert np.array_equal(out["total_extreme_count_annual"].sum("year").values, out["total_extreme_count"].values)
    assert np.allclose(out["total_extreme_precip_annual"].sum("year").values, out["total_extreme_precip"].values, rtol=1e-5)
    for k in TYPE_KEYS:
        assert np.allclose(out[f"{k}_precip_annual"].sum("year").values,
                           out[f"{k}_frac"].values * out["total_extreme_precip"].values, rtol=1e-4, atol=1e-5), k

    # a subset of the years: only those, in the order given; steps of the other years are left out of the per-year arrays
    sub = run(pr, thr, ds, annual_years=[2020])
    assert list(sub["year"].values) == [2020]
    assert np.array_equal(sub["total_extreme_precip_annual"].values, out["total_extreme_precip_annual"].isel(year=[1]).values)
    whole = [n for n in sub.data_vars if not n.endswith("_annual")]
    assert len(whole) == 2 + 2 * len(TYPE_KEYS)
    xr.testing.assert_identical(sub[whole], out[whole])


# ---------------------------------------------------------------------------------------------------------------------
def test_save_writes_year_int32_compressed_amounts_and_the_extra_attributes_last():
    pr, ds, thr = make_case(monthly_times())
    plain = run(pr, thr, ds)
    annual = run(pr, thr, ds, annual_years=YEARS)
    extra = {"threshold_file": "/some/thresholds.nc", "annual_years": "2019, 2020, 2021", "n_annual_years": 3,
             "min_year_coverage_days": 360}
    with tempfile.TemporaryDirectory() as d:
        plain_file, annual_file = f"{d}/plain.nc", f"{d}/annual.nc"
        sp.save_spatial_results(plain, plain_file, "testsrc", "2019-01-01", "2021-12-01", "P95")
        sp.save_spatial_results(annual, annual_file, "testsrc", "2019-01-01", "2021-12-01", "P95", extra_attrs=extra)

        with netCDF4.Dataset(annual_file) as nc, netCDF4.Dataset(plain_file) as nc0:
            assert nc.variables["year"].dtype == np.int32 and nc.variables["year"].dimensions == ("year",)
            assert list(nc.variables["year"][:]) == YEARS
            for name in ANNUAL_VARS:
                v = nc.variables[name]
                assert v.dimensions == ("year", "cell"), name
                assert v.filters()["zlib"] is True, name
                assert v.dtype == (np.int32 if name == "total_extreme_count_annual" else np.float32), name
            # the variables that were already there keep their storage exactly
            for name in nc0.variables:
                a, b = nc0.variables[name], nc.variables[name]
                assert (a.dtype, a.dimensions, a.filters()) == (b.dtype, b.dimensions, b.filters()), name
            # standard global attributes first and unchanged in order, the extra ones after them
            assert list(nc.ncattrs())[:len(nc0.ncattrs())] == list(nc0.ncattrs())
            assert list(nc.ncattrs())[len(nc0.ncattrs()):] == list(extra)
            assert nc.threshold_file == "/some/thresholds.nc" and nc.n_annual_years == 3
            assert nc.variables["total_extreme_precip_annual"].units == "mm/h"
            assert nc.variables["dc_precip_annual"].storm_type_code == 9

        with xr.open_dataset(plain_file, decode_cf=False) as a, xr.open_dataset(annual_file, decode_cf=False) as b:
            for name in a.variables:  # the stored values of every existing variable are equal
                assert np.array_equal(a[name].values, b[name].values, equal_nan=True), name


def test_default_threshold_file_follows_the_threshold_version_of_the_source():
    norm = os.path.normpath
    assert norm(sp.default_threshold_file({"source_name": "x"}, "/r/")) == "/r/extreme_precip/x_precip_percentiles_6h_hp8_v1.nc"
    assert (norm(sp.default_threshold_file({"source_name": "x", "threshold_version": "v9"}, "/r/"))
            == "/r/extreme_precip/x_precip_percentiles_6h_hp8_v9.nc")


def test_cli_defaults():
    old = sys.argv
    try:
        sys.argv = ["calc_stormtype_extreme_precip_spatial.py", "--catalog_source", "IR_IMERG"]
        args = sp.parse_args()
        sys.argv += ["--threshold_file", "/x/t.nc", "--no_annual", "--min_year_coverage_days", "300"]
        args2 = sp.parse_args()
    finally:
        sys.argv = old
    assert args.min_year_coverage_days == 360 and args.no_annual is False and args.threshold_file is None
    assert args2.min_year_coverage_days == 300 and args2.no_annual is True and args2.threshold_file == "/x/t.nc"


# ---------------------------------------------------------------------------------------------------------------------
def make_root(root, times, n_cells=48, seed=3):
    """A tiny pipeline data root: cof_masks/testsrc_cofmasks_hp8_v1.zarr (the masks and tot_pr on 48 HEALPix cells, nside 2),
    a default threshold file that makes nothing extreme and a config with the source. Returns (pr, ds, thresholds)."""
    pr, ds, thr = make_case(times, n_cells=n_cells, seed=seed)
    store = ds.assign(tot_pr=pr)
    for name in store.data_vars:
        store[name].attrs["grid_mapping"] = "crs"
    store = store.assign_coords(crs=xr.DataArray(np.zeros(0), dims=("crs",), attrs={
        "grid_mapping_name": "healpix", "healpix_nside": 2, "healpix_order": "nest"}))
    os.makedirs(f"{root}/cof_masks")
    os.makedirs(f"{root}/extreme_precip")
    store.to_zarr(f"{root}/cof_masks/testsrc_cofmasks_hp8_v1.zarr", consolidated=True)
    xr.Dataset({"pr_p90": ("cell", np.full(n_cells, np.inf, "float32"))}).to_netcdf(
        f"{root}/extreme_precip/testsrc_precip_percentiles_6h_hp8_v1.nc")
    (Path(root) / "config.yaml").write_text("TEST:\n  source_name: testsrc\nTESTV:\n  source_name: testsrc\n  threshold_version: v7\n")
    return pr, ds, thr


def run_main(root, output_dir, *extra_args, source="TEST"):
    old_argv, old_env = sys.argv, os.environ.get("COF_DATA_ROOT")
    sys.argv = ["calc_stormtype_extreme_precip_spatial.py", "--catalog_source", source, "--config_file", f"{root}/config.yaml",
                "--percentiles", "P90", "--output_dir", output_dir, "--n_workers", "2", "--batch_size", "7", *extra_args]
    os.environ["COF_DATA_ROOT"] = root
    try:
        sp.main()
    finally:
        sys.argv = old_argv
        if old_env is None:
            os.environ.pop("COF_DATA_ROOT", None)
        else:
            os.environ["COF_DATA_ROOT"] = old_env
    return f"{output_dir}/testsrc_stormtype_spatial_p90.nc"


def test_main_threshold_file_no_annual_and_min_year_coverage_days_end_to_end():
    times = monthly_times()
    with tempfile.TemporaryDirectory() as root:
        pr, ds, thr = make_root(root, times)
        alt = f"{root}/alt_thresholds.nc"
        xr.Dataset({"pr_p90": ("cell", thr)}).to_netcdf(alt)

        # the default threshold file (in the data root) makes nothing extreme
        f0 = run_main(root, f"{root}/out_default", "--min_year_coverage_days", "10")
        with xr.open_dataset(f0) as o:
            assert float(o["total_extreme_count"].sum()) == 0
            assert o.attrs["threshold_file"] == f"{root}/extreme_precip/testsrc_precip_percentiles_6h_hp8_v1.nc"

        # --threshold_file replaces it; 12 distinct days per year meet a 10-day gate, so the 3 years get per-year totals
        f1 = run_main(root, f"{root}/out_alt", "--threshold_file", alt, "--min_year_coverage_days", "10")
        with xr.open_dataset(f1) as o:
            assert o.attrs["threshold_file"] == alt
            assert o.attrs["annual_years"] == "2019, 2020, 2021" and o.attrs["n_annual_years"] == 3
            assert o.attrs["min_year_coverage_days"] == 10
            assert list(o["year"].values) == YEARS and ANNUAL_VARS <= set(o.data_vars)
            assert int(o["total_extreme_count"].sum()) > 0
            ref = run(pr, thr, ds, percentile_name="P90", annual_years=YEARS)
            assert np.array_equal(o["total_extreme_count"].values, ref["total_extreme_count"].values)
            assert np.allclose(o["total_extreme_precip_annual"].values, ref["total_extreme_precip_annual"].values)
            whole_record = {name: o[name].values for name in o.data_vars if not name.endswith("_annual")}

        # --no_annual: the same whole-record values, no per-year variables or attributes
        f2 = run_main(root, f"{root}/out_no_annual", "--threshold_file", alt, "--min_year_coverage_days", "10", "--no_annual")
        with xr.open_dataset(f2) as o:
            assert "year" not in o.coords and not any(n.endswith("_annual") for n in o.data_vars)
            assert "annual_years" not in o.attrs and o.attrs["threshold_file"] == alt
            for name, values in whole_record.items():
                assert np.array_equal(o[name].values, values, equal_nan=True), name

        # no --threshold_file: the default file follows the source's threshold_version (TESTV: v7), not the v1 file that
        # makes nothing extreme
        xr.Dataset({"pr_p90": ("cell", thr)}).to_netcdf(f"{root}/extreme_precip/testsrc_precip_percentiles_6h_hp8_v7.nc")
        f5 = run_main(root, f"{root}/out_version", "--min_year_coverage_days", "10", source="TESTV")
        with xr.open_dataset(f5) as o:
            assert o.attrs["threshold_file"] == f"{root}/extreme_precip/testsrc_precip_percentiles_6h_hp8_v7.nc"
            assert np.array_equal(o["total_extreme_count"].values, ref["total_extreme_count"].values)

        # a gate no year meets (12 days < 13): no per-year variables either, and no error
        f3 = run_main(root, f"{root}/out_strict", "--threshold_file", alt, "--min_year_coverage_days", "13")
        with xr.open_dataset(f3) as o:
            assert "year" not in o.coords and "annual_years" not in o.attrs
            for name, values in whole_record.items():
                assert np.array_equal(o[name].values, values, equal_nan=True), name


if __name__ == "__main__":
    test_qualifying_years_counts_distinct_days_and_needs_two_years()
    test_qualifying_years_with_a_cftime_calendar()
    test_whole_record_variables_are_identical_with_and_without_the_annual_option()
    test_whole_record_values_equal_a_plain_step_by_step_accumulation()
    test_annual_totals_match_direct_sums_and_close_the_whole_record()
    test_save_writes_year_int32_compressed_amounts_and_the_extra_attributes_last()
    test_default_threshold_file_follows_the_threshold_version_of_the_source()
    test_cli_defaults()
    test_main_threshold_file_no_annual_and_min_year_coverage_days_end_to_end()
    print("test_stormtype_annual: all checks passed")
