"""Tests of the per-calendar-year percentiles and interannual IQR added to scripts/calc_extreme_precip_thresholds.py (2026-09-21,
ported from calc_extreme_precip_thresholds_1h.py so the 6-hourly script can compute the same interannual variability for GSMaP
and IMERG):

  - calc_annual_precip_percentiles: returns (None, None) below the qualifying-year minimum; the coverage gate counts distinct
    calendar days, not time steps, so it is duration-agnostic; a year below --min_year_coverage_days is dropped; the annual
    values equal a per-year quantile computed directly with numpy; min_precip_threshold is applied the same way as
    calc_precip_percentiles's own filtering
  - calc_interannual_iqr: q25/q75/iqr equal numpy's quantile across the year dimension, for every cell
  - write_netcdf: with annual/iqr results, the file gains a year coordinate and the pr_annual_p*/pr_q25_p*/pr_q75_p*/pr_iqr_p*
    variables and the annual_years/n_annual_years attributes; without them (the --no_annual default off, or too few years), the
    file is exactly as before -- no new variables, no regression of the existing pr_p* values or attributes
  - calc_precip_percentiles (the existing all-record function) is untouched by the change: same values as a direct numpy quantile

Run:  python tests/test_extreme_precip_annual_iqr.py      (or pytest tests/)
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import calc_extreme_precip_thresholds as ct  # noqa: E402


def make_pr(n_cells=6, years=(2019, 2020, 2021), freq="6h", rng_seed=0, missing_year_days=None):
    """A synthetic (time, cell) precipitation DataArray: one 6-hourly time step every 6h across whole calendar years, values a
    reproducible mix of dry (< 0.1, to be filtered) and wet cells so the threshold and the quantiles are both exercised.
    missing_year_days: {year: days_to_keep} to make a year fall short of the coverage gate (keeps only the first N days)."""
    rng = np.random.default_rng(rng_seed)
    times = []
    for y in years:
        t = pd.date_range(f"{y}-01-01", f"{y}-12-31 18:00", freq=freq)
        if missing_year_days and y in missing_year_days:
            t = t[t.dayofyear <= missing_year_days[y]]
        times.append(t)
    time = pd.DatetimeIndex(np.concatenate(times))
    # cell 0: mostly dry (below the 0.1 mm/h threshold); the rest: gamma-distributed wet values, a different scale per cell so
    # the per-cell percentiles differ and per-year values differ (interannual spread to test the IQR against)
    data = np.stack([rng.gamma(shape=2.0, scale=0.05 + 0.3 * c / max(n_cells - 1, 1), size=len(time)).astype("float32")
                     for c in range(n_cells)], axis=1)
    data[:, 0] = rng.uniform(0, 0.05, size=len(time)).astype("float32")  # cell 0: always below threshold
    da = xr.DataArray(data, dims=["time", "cell"], coords={"time": time, "cell": np.arange(n_cells)})
    return da


def numpy_quantile_by_year(pr_da, year, percentiles, min_precip_threshold):
    sub = pr_da.sel(time=str(year)).values
    sub = np.where(sub >= min_precip_threshold, sub, np.nan)
    return {p: np.nanquantile(sub, p / 100.0, axis=0) for p in percentiles}


def test_annual_percentiles_match_numpy_and_respect_the_coverage_gate():
    pr = make_pr()
    results, years = ct.calc_annual_precip_percentiles(pr, percentiles=[90, 95], min_precip_threshold=0.1,
                                                        min_year_coverage_days=300, min_years=2)
    assert years == [2019, 2020, 2021]
    for p in (90, 95):
        assert results[p].dims == ("year", "cell") and results[p].sizes["year"] == 3
        for i, y in enumerate(years):
            want = numpy_quantile_by_year(pr, y, [p], 0.1)[p]
            got = results[p].isel(year=i).values
            assert np.allclose(got, want, equal_nan=True, atol=1e-5), (p, y)
    assert results[90].attrs["min_year_coverage_days"] == 300 and results[90].attrs["min_precip_threshold"] == "0.1 mm/h"

    # a year with only 60 days of data falls below the 300-day gate and is dropped, not just given fewer samples
    pr_short = make_pr(missing_year_days={2021: 60})
    results2, years2 = ct.calc_annual_precip_percentiles(pr_short, percentiles=[95], min_precip_threshold=0.1,
                                                         min_year_coverage_days=300, min_years=2)
    assert years2 == [2019, 2020]

    # below min_years, nothing is computed
    pr_one_year = make_pr(years=(2020,))
    r3, y3 = ct.calc_annual_precip_percentiles(pr_one_year, percentiles=[95], min_precip_threshold=0.1, min_years=2)
    assert r3 is None and y3 is None


def test_annual_percentiles_without_a_threshold_include_the_dry_cell():
    pr = make_pr()
    with_thresh, _ = ct.calc_annual_precip_percentiles(pr, percentiles=[95], min_precip_threshold=0.1, min_years=2)
    without_thresh, _ = ct.calc_annual_precip_percentiles(pr, percentiles=[95], min_precip_threshold=None, min_years=2)
    # cell 0 is uniform(0, 0.05): entirely excluded with the threshold (NaN), a small positive number without it
    assert np.isnan(with_thresh[95].isel(year=0, cell=0).values)
    assert 0 <= float(without_thresh[95].isel(year=0, cell=0).values) < 0.06


def test_interannual_iqr_matches_numpy_quantile_across_years():
    pr = make_pr()
    annual, years = ct.calc_annual_precip_percentiles(pr, percentiles=[90, 95], min_precip_threshold=0.1, min_years=2)
    iqr = ct.calc_interannual_iqr(annual, method="linear")
    for p in (90, 95):
        stacked = annual[p].values  # (year, cell)
        want_q25 = np.nanquantile(stacked, 0.25, axis=0)
        want_q75 = np.nanquantile(stacked, 0.75, axis=0)
        assert np.allclose(iqr[p]["q25"].values, want_q25, equal_nan=True, atol=1e-5)
        assert np.allclose(iqr[p]["q75"].values, want_q75, equal_nan=True, atol=1e-5)
        assert np.allclose(iqr[p]["iqr"].values, want_q75 - want_q25, equal_nan=True, atol=1e-5)
        assert iqr[p]["iqr"].attrs["n_years"] == 3 and iqr[p]["q25"].name == f"pr_q25_p{p}" and iqr[p]["iqr"].name == f"pr_iqr_p{p}"


def _fake_ds_p(n_cells):
    return xr.Dataset({}, coords={"cell": np.arange(n_cells), "lat": ("cell", np.linspace(-60, 60, n_cells)),
                                  "lon": ("cell", np.linspace(0, 300, n_cells)), "crs": 0})


def test_write_netcdf_without_annual_is_unchanged_from_before_the_port():
    pr = make_pr(n_cells=4)
    results = ct.calc_precip_percentiles(pr, percentiles=[90, 95], min_precip_threshold=0.1)
    results_computed = {p: d.compute() if hasattr(d, "compute") else d for p, d in results.items()}
    with tempfile.TemporaryDirectory() as d:
        path = f"{d}/out.nc"
        ct.write_netcdf(results_computed, _fake_ds_p(4), path, 8, "test_source", "2019-01-01", "2021-12-31", "6h", "linear",
                        min_precip_threshold=0.1)
        with xr.open_dataset(path) as ds:
            assert set(ds.data_vars) == {"pr_p90", "pr_p95"}
            assert "year" not in ds.coords and "annual_years" not in ds.attrs
            assert float(ds["pr_p95"].isel(cell=1).values) == float(results_computed[95].isel(cell=1).values)


def test_write_netcdf_with_annual_and_iqr_adds_the_expected_variables():
    pr = make_pr(n_cells=4)
    results = ct.calc_precip_percentiles(pr, percentiles=[95], min_precip_threshold=0.1)
    results_computed = {95: results[95].compute() if hasattr(results[95], "compute") else results[95]}
    annual, years = ct.calc_annual_precip_percentiles(pr, percentiles=[95], min_precip_threshold=0.1, min_years=2)
    iqr = ct.calc_interannual_iqr(annual)
    with tempfile.TemporaryDirectory() as d:
        path = f"{d}/out.nc"
        ct.write_netcdf(results_computed, _fake_ds_p(4), path, 8, "test_source", "2019-01-01", "2021-12-31", "6h", "linear",
                        min_precip_threshold=0.1, annual_results=annual, iqr_results=iqr, years=years,
                        min_year_coverage_days=300)
        with xr.open_dataset(path) as ds:
            assert set(ds.data_vars) == {"pr_p95", "pr_annual_p95", "pr_q25_p95", "pr_q75_p95", "pr_iqr_p95"}
            assert list(ds["year"].values) == years and ds.sizes["year"] == 3
            assert ds.attrs["annual_years"] == "2019, 2020, 2021" and ds.attrs["n_annual_years"] == 3
            assert ds.attrs["min_year_coverage_days"] == 300
            assert ds["pr_annual_p95"].dims == ("year", "cell")
            assert np.allclose(ds["pr_iqr_p95"].values, (ds["pr_q75_p95"] - ds["pr_q25_p95"]).values, equal_nan=True)
            assert ds["pr_annual_p95"].attrs["cell_methods"] == "time: quantile (interval: 1 year)"


def test_calc_precip_percentiles_is_unaffected_by_the_port():
    """The existing all-record function must still match a direct numpy quantile, unchanged by the functions added next to it."""
    pr = make_pr(n_cells=5)
    results = ct.calc_precip_percentiles(pr, percentiles=[90, 95], min_precip_threshold=0.1)
    masked = np.where(pr.values >= 0.1, pr.values, np.nan)
    for p in (90, 95):
        want = np.nanquantile(masked, p / 100.0, axis=0)
        got = results[p].values if not hasattr(results[p], "compute") else results[p].compute().values
        assert np.allclose(got, want, equal_nan=True, atol=1e-5)


def test_cli_defaults_and_gsmap_config_entry():
    out = subprocess.run([sys.executable, str(REPO / "scripts" / "calc_extreme_precip_thresholds.py"), "--help"],
                         capture_output=True, text=True).stdout
    text = " ".join(out.split())
    assert "--no_annual" in text and "--min_year_coverage_days" in text and "default: 300" in text

    import yaml
    cfg = yaml.safe_load(open(REPO / "config" / "config_sources.yaml"))
    assert "GSMAP" in cfg
    g = cfg["GSMAP"]
    assert g["source_name"] == "GSMaPv8" and g["varname_precip_liq"] == "hourlyPrecipRate" and g["pr_convert_factor"] == 1.0
    assert g.get("varname_precip_ice") is None


if __name__ == "__main__":
    test_annual_percentiles_match_numpy_and_respect_the_coverage_gate()
    test_annual_percentiles_without_a_threshold_include_the_dry_cell()
    test_interannual_iqr_matches_numpy_quantile_across_years()
    test_write_netcdf_without_annual_is_unchanged_from_before_the_port()
    test_write_netcdf_with_annual_and_iqr_adds_the_expected_variables()
    test_calc_precip_percentiles_is_unaffected_by_the_port()
    test_cli_defaults_and_gsmap_config_entry()
    print("test_extreme_precip_annual_iqr: all checks passed")
