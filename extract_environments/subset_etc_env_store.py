#!/usr/bin/env python
"""
Build the environment store of a NEW ETC track file from the store of the OLD track file, without extracting from the catalog again.

Use it when a re-tracking only removes storm points (or renumbers storms): the value of an environment variable around a storm point depends only on the
point's time and position (and the catalog), not on the storm ID or on the other points, so the store of the new points is a selection of rows of the old
store. The rows are matched by (time, lon, lat) of the points; every new point must be in the old list (else the script fails and the variable has to be
extracted again, --extract-env of run_etc_pipeline.py). The new store is written with the extractor's own save_to_zarr(), with the storm IDs, grid IDs,
positions and times of the NEW track file, so that it is the store that extract_etc_2d_vars.py would have written (same layout, chunks and compressor).

A point that two storms share has one row in the old store per storm; the first is used (the values are the same, they depend on time and position only).
The point list is the one of the extraction: rows with a position, |lat| <= 90 - radius, in the order of the track file.

Checks before writing: the old store belongs to the old track file (same number of points, and the same storm IDs, grid IDs, positions and times row by
row), every new point is found, and the matched rows have the grid ID and position of the new point. After writing the store is read back and compared.

Usage:
  python subset_etc_env_store.py --src-store OLD/etc_2d_tas_all_all.zarr --dst-store NEW/etc_2d_tas_all_all.zarr \\
      --old-track-file OLD_TRACKS.txt --new-track-file NEW_TRACKS.txt [--radius 20]

Structured meshes (ERA5, lon_id/lat_id) are not supported: extract those again.
Exit status: 0 written; 1 refused or failed (nothing is left at --dst-store).

Author: Zhe Feng | zhe.feng@pnnl.gov
"""
import argparse
import contextlib
import io
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
for p in (REPO, Path(__file__).resolve().parent):
    sys.path.insert(0, str(p))
from src.env_extract_utilities import parse_etc_track_file          # noqa: E402

STALE_ATTRS = ("n_points", "n_points_time_missing", "n_points_failed", "extraction_info")


def point_list(track_file, radius):
    """Storm points as the extraction sees them: with a position, |lat| <= 90 - radius, in file order (index reset)."""
    with contextlib.redirect_stdout(io.StringIO()):
        df = parse_etc_track_file(track_file, unstructured_mesh=True)
    df = df.dropna(subset=["lat", "lon"])
    df = df[df["lat"].between(-90 + radius, 90 - radius)].reset_index(drop=True)
    return df


def point_keys(df):
    """(time in ns, lon, lat) of every point, positions rounded to 1e-6 degrees (the text file has 6 digits)."""
    t = pd.to_datetime(df["base_time"]).values.astype("datetime64[ns]").astype("int64")
    return list(zip(t.tolist(), np.round(df["lon"].values, 6).tolist(), np.round(df["lat"].values, 6).tolist()))


def match_points(old_df, new_df):
    """Row of the old list for every new point: index array. Raises ValueError when a new point is not in the old list or its grid ID differs."""
    first = {}
    for i, k in enumerate(point_keys(old_df)):
        first.setdefault(k, i)
    pos = np.array([first.get(k, -1) for k in point_keys(new_df)], dtype=np.int64)
    n_missing = int((pos < 0).sum())
    if n_missing:
        raise ValueError(f"{n_missing} of {len(new_df)} new storm points are not in the old track file (time, lon, lat): extract the variable again")
    if not (old_df["grid_id"].values[pos] == new_df["grid_id"].values).all():
        raise ValueError("the grid IDs of matching points differ between the old and the new track file")
    return pos


def check_store_belongs_to(src, old_df):
    """The old store must be the extraction of the old point list: same count and the same metadata row by row."""
    if src.sizes["time"] != len(old_df):
        raise ValueError(f"the old store has {src.sizes['time']} storm points, the old track file gives {len(old_df)}")
    if not (src["storm_id"].values == old_df["storm_id"].values).all():
        raise ValueError("the storm IDs of the old store differ from the old track file")
    if not (src["grid_id"].values == old_df["grid_id"].values).all():
        raise ValueError("the grid IDs of the old store differ from the old track file")
    if not (np.allclose(src["storm_lat"].values, old_df["lat"].values) and np.allclose(src["storm_lon"].values, old_df["lon"].values)):
        raise ValueError("the positions of the old store differ from the old track file")
    t_store = pd.to_datetime(src["time"].values).values.astype("datetime64[ns]")
    if not (t_store == pd.to_datetime(old_df["base_time"]).values.astype("datetime64[ns]")).all():
        raise ValueError("the times of the old store differ from the old track file")


def subset_store(src_store, dst_store, old_track_file, new_track_file, radius=None):
    import xarray as xr
    import zarr
    from extract_etc_2d_vars import save_to_zarr

    t0 = time.time()
    src = xr.open_zarr(src_store)
    if "lon_id" in src.data_vars:
        raise ValueError("structured mesh stores (lon_id/lat_id) are not supported: extract the variable again")
    radius = float(src.attrs["radius"]) if radius is None else float(radius)
    if abs(radius - float(src.attrs["radius"])) > 1e-9:
        raise ValueError(f"radius {radius} differs from the store's {src.attrs['radius']}")
    variables = [v for v in src.data_vars if src[v].dims == ("time", "y", "x")]
    if len(variables) != 1:
        raise ValueError(f"expected one (time, y, x) variable in {src_store}, found {variables}")
    var = variables[0]

    old_df, new_df = point_list(old_track_file, radius), point_list(new_track_file, radius)
    check_store_belongs_to(src, old_df)
    pos = match_points(old_df, new_df)
    print(f"{os.path.basename(src_store)}: {len(new_df)} of {len(old_df)} old points kept ({len(old_df) - len(set(pos.tolist()))} old rows not used)", flush=True)

    data = src[var].isel(time=pos).values                                # (n_new, y, x), float32 as extracted
    attrs = {k: v for k, v in src[var].attrs.items() if k not in STALE_ATTRS}
    chunk_size = int(src[var].encoding.get("chunks", (1000,))[0])
    out = str(dst_store)[:-len(".zarr")] if str(dst_store).endswith(".zarr") else str(dst_store)
    try:
        save_to_zarr(data, np.array(pd.to_datetime(new_df["base_time"])), new_df["storm_id"].values, new_df["grid_id"].values,
                     new_df["lat"].values, new_df["lon"].values, src["x"].values, src["y"].values, var, out,
                     radius, float(src.attrs["lon_res"]), float(src.attrs["lat_res"]), chunk_size=chunk_size, unstructured_mesh=True, var_attrs=attrs)
        zarr_path = out + ".zarr"
        g = zarr.open_group(zarr_path, mode="r+")
        g.attrs["reused_from"] = str(src_store)
        g.attrs["reused_old_track_file"] = str(old_track_file)
        g.attrs["reused_new_track_file"] = str(new_track_file)
        g.attrs["reused_n_points"] = int(len(new_df))
        zarr.consolidate_metadata(zarr_path)
        # read back
        chk = xr.open_zarr(zarr_path)
        assert chk.sizes["time"] == len(new_df), "wrong number of points written"
        assert (chk["storm_id"].values == new_df["storm_id"].values).all() and (chk["grid_id"].values == new_df["grid_id"].values).all()
        rows = np.unique(np.linspace(0, len(new_df) - 1, 7).astype(int))
        a, b = chk[var].isel(time=rows).values, src[var].isel(time=pos[rows]).values
        assert ((a == b) | (np.isnan(a) & np.isnan(b))).all(), "the values read back differ from the selected rows"
    except BaseException:
        shutil.rmtree(out + ".zarr", ignore_errors=True)
        raise
    print(f"  wrote {out}.zarr ({var}, {len(new_df)} points) in {time.time() - t0:.0f} s", flush=True)
    return out + ".zarr"


def main():
    ap = argparse.ArgumentParser(description="Select the rows of an old ETC environment store that belong to the points of a new track file.")
    ap.add_argument("--src-store", required=True, help="old store, etc_2d_<var>_all_all.zarr")
    ap.add_argument("--dst-store", required=True, help="new store to write (must not exist)")
    ap.add_argument("--old-track-file", required=True, help="the track file the old store was extracted for")
    ap.add_argument("--new-track-file", required=True)
    ap.add_argument("--radius", type=float, default=None, help="extraction radius in degrees (default: the store's)")
    args = ap.parse_args()
    if os.path.exists(args.dst_store):
        print(f"ERROR: {args.dst_store} exists; nothing is overwritten", flush=True)
        return 1
    try:
        os.makedirs(os.path.dirname(os.path.abspath(args.dst_store)), exist_ok=True)
        subset_store(args.src_store, args.dst_store, args.old_track_file, args.new_track_file, args.radius)
    except Exception as exc:                                              # noqa: BLE001
        print(f"ERROR: {type(exc).__name__}: {exc}", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
