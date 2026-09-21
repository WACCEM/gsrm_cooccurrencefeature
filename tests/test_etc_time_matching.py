"""Tests of the exact time matching in extract_environments/extract_etc_2d_vars.py (extract_etc_2d_variable).

  - a track time that exists in the source gives that frame (checked against a value that encodes the frame and the cell)
  - a track time without a frame gives a NaN slab, is counted, keeps its metadata, and is NEVER replaced by the nearest frame
    (a time after the end or before the start of the record used to return the last / first frame)
  - the chunk cache (whole time chunks read once) gives exactly the same slabs as reading frame by frame
  - a repeated time stamp in the source uses the first frame
  - cof_source_name maps the catalog model names to the COF store names

Run:  python tests/test_etc_time_matching.py      (or pytest tests/)
"""
import contextlib
import io
import sys
from pathlib import Path

import healpy as hp
import numpy as np
import pandas as pd
import xarray as xr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "extract_environments"))
sys.path.insert(0, str(REPO))
import extract_etc_2d_vars as ex  # noqa: E402

NSIDE, RADIUS, RES = 8, 2.0, 1.0
T0 = np.datetime64("2020-01-01T00:00", "ns")
TIMES = T0 + np.arange(10) * np.timedelta64(6, "h")


def source(chunk_time=None, times=TIMES):
    n = len(times)
    values = (1000.0 * np.arange(n)[:, None] + np.arange(hp.nside2npix(NSIDE))[None, :]).astype("float32")
    da = xr.DataArray(values, dims=("time", "cell"), coords={"time": times, "cell": np.arange(values.shape[1])}, attrs={"units": "mm h-1"})
    return da.chunk({"time": chunk_time, "cell": -1}) if chunk_time else da


def storms():
    rows = [  # storm_id, lon, lat, time
        (1, 10.0, 20.0, TIMES[2]), (2, 100.0, -30.0, TIMES[2]), (3, 45.0, 10.0, TIMES[7]),
        (4, 0.0, 0.0, TIMES[3] + np.timedelta64(3, "h")),       # between two frames
        (5, 30.0, 30.0, TIMES[-1] + np.timedelta64(6, "h")),    # after the end of the record
        (6, 60.0, 40.0, TIMES[0] - np.timedelta64(6, "h")),     # before the start
    ]
    return pd.DataFrame([{"storm_id": i, "lon": lo, "lat": la, "base_time": t, "grid_id": 100 + i} for i, lo, la, t in rows])


def expected(frame, lon, lat):
    lons, lats = np.meshgrid(np.arange(lon - RADIUS, lon + RADIUS + RES, RES), np.arange(lat - RADIUS, lat + RADIUS + RES, RES))
    return (1000.0 * frame + hp.ang2pix(NSIDE, lons, lats, nest=True, lonlat=True)).astype("float32")


def run(variable_data, df=None):
    df = storms() if df is None else df
    with contextlib.redirect_stdout(io.StringIO()):
        return df, ex.extract_etc_2d_variable(df, variable_data, xr.Dataset(), NSIDE, radius=RADIUS, lon_res=RES, lat_res=RES, variable_name="v")


def check(df, res):
    out, time_arr, ids, gids, lats, lons, xc, yc, attrs = res
    frame_of = {1: 2, 2: 2, 3: 7}
    for k, row in df.iterrows():
        if row["storm_id"] in frame_of:
            np.testing.assert_array_equal(out[k], expected(frame_of[row["storm_id"]], row["lon"], row["lat"]))
        else:
            assert np.isnan(out[k]).all(), f"storm {row['storm_id']} has no frame at its time: the slab must be NaN, never the nearest frame"
        # the metadata is kept for every point, also for the NaN ones
        assert ids[k] == row["storm_id"] and gids[k] == row["grid_id"] and pd.Timestamp(time_arr[k]) == pd.Timestamp(row["base_time"])
    assert attrs["n_points"] == 6 and attrs["n_points_time_missing"] == 3 and attrs["n_points_failed"] == 0, attrs
    assert attrs["units"] == "mm h-1" and "exact" in attrs["time_match"]
    return out


def test_exact_matching_and_cache():
    df, r_numpy = run(source(chunk_time=None))                 # in memory: no chunks, no cache
    out_numpy = check(df, r_numpy)
    df, r_cache = run(source(chunk_time=4))                    # dask, time chunks 4,4,2: the cache is used
    out_cache = check(df, r_cache)
    saved = ex.MAX_CACHED_CHUNK_BYTES
    ex.MAX_CACHED_CHUNK_BYTES = 1                              # dask but too big to cache: frame by frame
    try:
        df, r_frame = run(source(chunk_time=4))
        out_frame = check(df, r_frame)
    finally:
        ex.MAX_CACHED_CHUNK_BYTES = saved
    np.testing.assert_array_equal(out_numpy, out_cache)
    np.testing.assert_array_equal(out_numpy, out_frame)


def test_repeated_time_stamp_uses_first_frame():
    times = TIMES.copy(); times[5] = times[4]                  # frames 4 and 5 carry the same stamp
    df = pd.DataFrame([{"storm_id": 9, "lon": 10.0, "lat": 20.0, "base_time": TIMES[4], "grid_id": 109}])
    _, res = run(source(times=times), df)
    np.testing.assert_array_equal(res[0][0], expected(4, 10.0, 20.0))       # the first of the two frames


def test_cof_source_name():
    f = ex.cof_source_name
    assert f("scream2D_hrly") == "scream" and f("scream_ne120_inst") == "scream" and f("scream_ne120") == "scream"
    assert f("era5_3h") == "IMERGv7"
    assert f("nicam_gl11_shifted") == "nicam_gl11" and f("nicam_gl11") == "nicam_gl11"
    assert f("icon_d3hp003") == "icon_d3hp003"
    assert f("um_glm_n2560_RAL3p3") == "um_glm_n2560_RAL3p3"
    assert f("casesm2_10km_nocumulus") == "casesm2_10km_nocumulus"


if __name__ == "__main__":
    test_exact_matching_and_cache(); test_repeated_time_stamp_uses_first_frame(); test_cof_source_name()
    print("test_etc_time_matching: all checks passed")
