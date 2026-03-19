#!/usr/bin/env python3
"""
check_etc_storm.py — Look up an ETC storm by ID in the stitched_nodes text file.

The track file format (ERA5 structured, others unstructured HEALPix):
  start  NUM_TIMESTEPS  YEAR  MONTH  DAY  HOUR
    LON_ID  LAT_ID  LON  LAT  ...extra cols...  YEAR  MONTH  DAY  HOUR  (structured)
    GRID_ID  LON  LAT  ...extra cols...  YEAR  MONTH  DAY  HOUR            (unstructured)

Storm IDs are sequential: the Nth "start" line → storm_id = N.

Usage examples:
  python check_etc_storm.py --storm_id 1 --trackfile era5.etc_stitched_nodes.txt
  python check_etc_storm.py --storm_id 535 --unstructured
  python check_etc_storm.py --storm_id 535 --check_cof_mask
  python check_etc_storm.py --list_range 530 540
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd

# Add parent directory to path to import from src
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.env_extract_utilities import parse_etc_track_file


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# COF mask check
# ---------------------------------------------------------------------------

def check_cof_mask(storm_df, cof_zarr_path):
    """
    For each timestep of the storm, print the ETC IDs visible in the COF mask
    within a small lat/lon box around the storm centre.

    Parameters
    ----------
    storm_df  : DataFrame rows for one storm (must have lat, lon, base_time)
    cof_zarr_path : path to e.g. IMERGv7_cofmasks_hp8_v1.zarr
    """
    try:
        import zarr, healpy as hp
        import xarray as xr
    except ImportError as e:
        print(f"  [skip] COF mask check requires zarr/healpy/xarray: {e}")
        return

    print(f"\n{'='*60}")
    print(f"COF mask check: {cof_zarr_path}")
    print(f"{'='*60}")

    z = zarr.open(cof_zarr_path, "r")

    # Decode times
    time_raw = z["time"][:]
    time_attrs = dict(z["time"].attrs)
    units = time_attrs.get("units", "hours since 2019-01-01")
    cal   = time_attrs.get("calendar", "proleptic_gregorian")
    try:
        import cftime
        times_cf = cftime.num2date(time_raw, units=units, calendar=cal)
        times_pd = pd.DatetimeIndex([pd.Timestamp(t.year, t.month, t.day, t.hour) for t in times_cf])
    except Exception:
        ref = pd.Timestamp(units.split("since")[-1].strip())
        times_pd = ref + pd.to_timedelta(time_raw, unit="h")

    # HEALPix order — try to infer nside from array length
    npix = z["etc_mask"].shape[1]
    nside = hp.npix2nside(npix)
    nest  = True   # assume nested; adjust if ring

    print(f"  HEALPix nside={nside}, npix={npix}\n")

    for _, row in storm_df.iterrows():
        t = row["base_time"]
        tidx = np.searchsorted(times_pd, t)
        if tidx >= len(times_pd) or times_pd[tidx] != t:
            continue

        etc_slice = z["etc_mask"][tidx, :]

        # Find pixels within ~5° of storm centre
        theta = np.radians(90 - row["lat"])
        phi   = np.radians(row["lon"])
        vec   = hp.ang2vec(theta, phi)
        radius_rad = np.radians(5.0)
        pix_near = hp.query_disc(nside, vec, radius_rad, nest=nest)

        ids_near = np.unique(etc_slice[pix_near])
        ids_near = ids_near[ids_near > 0].astype(int)

        print(f"  {t}  lat={row['lat']:.2f} lon={row['lon']:.2f}  "
              f"COF etc IDs within 5°: {ids_near.tolist()}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    default_trackfile = "/pscratch/sd/w/wcmca1/hackathon/etc_tracks/era5.etc_stitched_nodes.txt"
    default_cof_zarr  = "/pscratch/sd/w/wcmca1/hackathon/cof_masks/IMERGv7_cofmasks_hp8_v1.zarr"

    parser = argparse.ArgumentParser(
        description="Look up an ETC storm by sequential ID in a stitched_nodes text file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--trackfile", default=default_trackfile,
                        help="Path to etc_stitched_nodes.txt  (default: ERA5)")
    parser.add_argument("--storm_id", type=int, default=None,
                        help="Sequential storm ID to look up (1-based, Nth 'start' line)")
    parser.add_argument("--list_range", nargs=2, type=int, metavar=("START", "END"),
                        help="Print summary for all storm IDs in [START, END] range")
    parser.add_argument("--unstructured", action="store_true",
                        help="Use unstructured (HEALPix) column layout. "
                             "Default: structured lat/lon (ERA5).")
    parser.add_argument("--check_cof_mask", action="store_true",
                        help="Also query the IMERG COF mask zarr for nearby ETC IDs "
                             "(useful to cross-check sequential IDs vs mask pixel IDs)")
    parser.add_argument("--cof_zarr", default=default_cof_zarr,
                        help="Path to COF mask zarr used with --check_cof_mask")
    args = parser.parse_args()

    if args.storm_id is None and args.list_range is None:
        parser.error("Provide either --storm_id or --list_range")

    # ---------- parse file ----------
    print(f"Reading: {args.trackfile}")
    df = parse_etc_track_file(args.trackfile, unstructured_mesh=args.unstructured)
    total = df["storm_id"].max()
    print(f"Total storms (start lines): {total}")
    print(f"Time range: {df['base_time'].min()} → {df['base_time'].max()}\n")

    # ---------- handle --list_range ----------
    if args.list_range:
        s0, s1 = args.list_range
        sub = df[df["storm_id"].between(s0, s1)]
        if sub.empty:
            print(f"No storms in range {s0}–{s1}")
            return
        print(f"{'ID':>6}  {'npts':>5}  {'start':20s}  {'end':20s}  "
              f"{'lat0':>7}  {'lon0':>7}")
        print("-" * 70)
        for sid, grp in sub.groupby("storm_id"):
            grp_sorted = grp.sort_values("base_time")
            print(f"{sid:6d}  {len(grp):5d}  "
                  f"{str(grp_sorted['base_time'].iloc[0])[:19]:20s}  "
                  f"{str(grp_sorted['base_time'].iloc[-1])[:19]:20s}  "
                  f"{grp_sorted['lat'].iloc[0]:7.2f}  "
                  f"{grp_sorted['lon'].iloc[0]:7.2f}")
        return

    # ---------- handle --storm_id ----------
    sid = args.storm_id
    storm_df = df[df["storm_id"] == sid].sort_values("base_time").reset_index(drop=True)

    if storm_df.empty:
        print(f"Storm ID {sid} not found. Valid range: 1–{total}")
        sys.exit(1)

    print(f"{'='*60}")
    print(f"Storm ID : {sid}  (sequential, Nth 'start' line)")
    print(f"Points   : {len(storm_df)}")
    nt = storm_df["num_timesteps"].iloc[0]
    print(f"Declared num_timesteps in file: {nt}")
    print(f"{'='*60}")

    # column layout
    if args.unstructured:
        header = f"{'#':>4}  {'base_time':20s}  {'lat':>8}  {'lon':>9}  {'grid_id':>10}"
    else:
        header = f"{'#':>4}  {'base_time':20s}  {'lat':>8}  {'lon':>9}  {'lon_id':>7}  {'lat_id':>7}"
    print(header)
    print("-" * len(header))

    for i, row in storm_df.iterrows():
        if args.unstructured:
            print(f"{i+1:4d}  {str(row['base_time'])[:19]:20s}  "
                  f"{row['lat']:8.3f}  {row['lon']:9.3f}  "
                  f"{int(row['grid_id']):10d}")
        else:
            print(f"{i+1:4d}  {str(row['base_time'])[:19]:20s}  "
                  f"{row['lat']:8.3f}  {row['lon']:9.3f}  "
                  f"{int(row['lon_id']):7d}  {int(row['lat_id']):7d}")

    if args.check_cof_mask:
        check_cof_mask(storm_df, args.cof_zarr)


if __name__ == "__main__":
    main()
