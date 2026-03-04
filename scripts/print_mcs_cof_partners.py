#!/usr/bin/env python3
"""
Print the AR and ETC partner track IDs for a given MCS track number.

Usage
-----
python print_mcs_cof_partners.py --source scream --mcs-id 476
python print_mcs_cof_partners.py --source scream --mcs-id 1024
python print_mcs_cof_partners.py --nc /path/to/custom_2d.nc --mcs-id 476

MCS ID convention
-----------------
Use the COF mask 1-based track ID (i.e., the value stored in the COF zarr
mask arrays).  Internally this maps to stats track index = mcs_id - 1.
"""

import argparse
import os
from datetime import datetime, timezone

import numpy as np
import xarray as xr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Print COF partner AR/ETC track IDs for one MCS track."
    )
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--source', default='scream',
                       help='Source name used to locate the default 2D netCDF '
                            '(e.g. scream, icon_d3hp003).  Ignored when --nc is given.')
    group.add_argument('--nc', dest='nc_path',
                       help='Direct path to the *_mcs_cof_tracks_2d.nc file.')
    parser.add_argument('--mcs-id', type=int, required=True,
                        help='MCS COF mask track ID (1-based, as stored in the '
                             'COF zarr mask arrays).')
    parser.add_argument('--cof-dir', default='/pscratch/sd/w/wcmca1/hackathon/cof_masks',
                        help='Root directory for COF output '
                             '(default: /pscratch/sd/w/wcmca1/hackathon/cof_masks).')
    return parser.parse_args()


def main():
    args = parse_args()

    if args.nc_path:
        nc_path = args.nc_path
    else:
        nc_path = os.path.join(
            args.cof_dir, 'stats', f'{args.source}_mcs_cof_tracks_2d.nc'
        )

    if not os.path.exists(nc_path):
        raise FileNotFoundError(f"2D tracks netCDF not found: {nc_path}")

    mcs_cof_id  = args.mcs_id
    track_idx   = mcs_cof_id - 1   # 0-based stats index

    ds = xr.open_dataset(nc_path, mask_and_scale=False, decode_times=False)
    n_tracks = ds.sizes['tracks']

    if track_idx < 0 or track_idx >= n_tracks:
        raise ValueError(
            f"MCS COF ID {mcs_cof_id} (index {track_idx}) is out of range "
            f"(file has {n_tracks} tracks)."
        )

    base_time = ds['base_time'].values[track_idx].astype(float)
    ar_ids    = ds['ar_tracknum'].values[track_idx]
    etc_ids   = ds['etc_tracknum'].values[track_idx]

    # Optional 3-way arrays (may not exist in older outputs)
    ar_3way  = ds['ar_tracknum_3way'].values[track_idx]  if 'ar_tracknum_3way'  in ds else None
    etc_3way = ds['etc_tracknum_3way'].values[track_idx] if 'etc_tracknum_3way' in ds else None

    valid = base_time > 0

    print(f"\nMCS COF track ID : {mcs_cof_id}  (stats index {track_idx})")
    print(f"File             : {nc_path}")
    print(f"Total time steps : {len(base_time)}")
    print(f"Valid time steps : {valid.sum()}")

    has_3way = (ar_3way is not None) and (etc_3way is not None)
    if has_3way:
        header = f"{'Step':>5}  {'UTC time (base_time)':22}  {'AR':>6}  {'ETC':>6}  {'AR(3way)':>10}  {'ETC(3way)':>10}"
        sep    = "-" * 70
    else:
        header = f"{'Step':>5}  {'UTC time (base_time)':22}  {'AR':>6}  {'ETC':>6}"
        sep    = "-" * 46

    print()
    print(header)
    print(sep)

    for i in np.where(valid)[0]:
        t       = datetime.fromtimestamp(float(base_time[i]), tz=timezone.utc).strftime('%Y-%m-%d %H:%M')
        ar_str  = str(int(ar_ids[i]))  if int(ar_ids[i])  > 0 else '-'
        etc_str = str(int(etc_ids[i])) if int(etc_ids[i]) > 0 else '-'

        if has_3way:
            ar3_str  = str(int(ar_3way[i]))  if int(ar_3way[i])  > 0 else '-'
            etc3_str = str(int(etc_3way[i])) if int(etc_3way[i]) > 0 else '-'
            print(f"{i:>5}  {t:22}  {ar_str:>6}  {etc_str:>6}  {ar3_str:>10}  {etc3_str:>10}")
        else:
            print(f"{i:>5}  {t:22}  {ar_str:>6}  {etc_str:>6}")

    ds.close()
    print()


if __name__ == '__main__':
    main()
