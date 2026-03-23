"""
Combine MCS track statistics NetCDF files across multiple years into a single NetCDF output.

- Shifts the tracks coordinate in each year so that MCS track IDs are unique and sequential across years.
- Preserves the original data structure and missing values.
- Uses minimal compression for the output file.

Usage:
  python combine_mcs_track_stats.py [--overwrite]

Options:
  --overwrite   Overwrite the output NetCDF file if it already exists.

Output:
  Writes a single NetCDF file combining all input years to:
    /pscratch/sd/w/wcmca1/hackathon/mcs/IMERGv7_2019_2021/stats/mcs_tracks_final_startdate_enddate.nc
"""

import os
import glob
import re
import xarray as xr
import numpy as np
import argparse
import time

def main():
    parser = argparse.ArgumentParser(description="Combine MCS track stats NetCDF files across years.")
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output file if it exists')
    args = parser.parse_args()

    stats_dir = "/pscratch/sd/w/wcmca1/hackathon/mcs/IMERGv7_2019_2021/stats"
    pattern = os.path.join(stats_dir, "mcs_tracks_final_20??????.????_20??????.????.nc")

    files = sorted(glob.glob(pattern))
    if not files:
        raise RuntimeError("No input files found!")

    print(f"Found {len(files)} files:")
    for f in files:
        print(f"  {f}")

    # Extract start/end dates for output filename
    date_re = re.compile(r"mcs_tracks_final_(\d{8}\.\d{4})_(\d{8}\.\d{4})\.nc")
    all_starts, all_ends = [], []
    for f in files:
        m = date_re.search(os.path.basename(f))
        if m:
            all_starts.append(m.group(1))
            all_ends.append(m.group(2))
    startdate = min(all_starts)
    enddate = max(all_ends)
    out_path = os.path.join(stats_dir, f"mcs_tracks_final_{startdate}_{enddate}.nc")

    if os.path.exists(out_path) and not args.overwrite:
        print(f"Output file already exists: {out_path}")
        print("Use --overwrite to overwrite the output file.")
        return

    t0 = time.time()
    datasets = []
    track_offset = 0
    for i, f in enumerate(files):
        t1 = time.time()
        print(f"Reading file {i+1}/{len(files)}: {f}")
        ds = xr.open_dataset(f, mask_and_scale=False, decode_times=False)
        ntracks = ds.sizes["tracks"]
        ds = ds.assign_coords(tracks=ds.tracks + track_offset)
        datasets.append(ds)
        print(f"  Tracks in file: {ntracks}, track_offset now: {track_offset}")
        track_offset += ntracks
        print(f"  Done reading file {i+1} ({time.time() - t1:.2f} s elapsed)")

    print(f"Concatenating {len(datasets)} datasets...")
    t2 = time.time()
    combined = xr.concat(datasets, dim="tracks")
    print(f"Concatenation done ({time.time() - t2:.2f} s)")

    print(f"Writing output to {out_path} ...")
    t3 = time.time()
    encoding = {v: {"zlib": True, "complevel": 1} for v in combined.data_vars}
    combined.to_netcdf(out_path, encoding=encoding)
    print(f"Combined file written: {out_path}")
    print(f"Total time: {time.time() - t0:.2f} s")

if __name__ == "__main__":
    main()
