"""
Combine MCS pixel-level mask Zarr files across multiple years into a single Zarr output.

- Shifts the mcs_mask values in each year so that MCS IDs are unique and sequential across years.
- Preserves the original chunking and compression from the input files.
- Uses Dask for parallel writing to efficiently handle large datasets.
- Only mcs_mask is shifted; ccs_mask is concatenated as-is.

Usage:
  python combine_mcs_mask_zarr.py [--overwrite]

Options:
  --overwrite   Overwrite the output Zarr file if it already exists.

Output:
  Writes a single Zarr file combining all input years to:
    /pscratch/sd/w/wcmca1/hackathon/mcs/IMERGv7_2019_2021/mcstracking/IMERGv7_hrly_mcsmask_hp8_v1.zarr
"""

import os
import glob
import re
import xarray as xr
import numpy as np
import argparse
import time
import dask
from dask.diagnostics import ProgressBar

def main():
    parser = argparse.ArgumentParser(description="Combine MCS mask Zarr files across years.")
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output file if it exists')
    args = parser.parse_args()

    # Input/output paths
    mask_dir = "/pscratch/sd/w/wcmca1/hackathon/mcs/IMERGv7_2019_2021/mcstracking"
    pattern = os.path.join(mask_dir, "IMERGv7_hrly_mcsmask_hp8_v1_20??????.????_20??????.????.zarr")
    out_path = os.path.join(mask_dir, "IMERGv7_hrly_mcsmask_hp8_v1.zarr")

    # Find and sort files by start date
    files = sorted(glob.glob(pattern))
    if not files:
        raise RuntimeError("No input files found!")
    
    print(f"Found {len(files)} files:")
    for f in files:
        print(f"  {f}")

    if os.path.exists(out_path) and not args.overwrite:
        print(f"Output file already exists: {out_path}")
        print("Use --overwrite to overwrite the output file.")
        return

    t0 = time.time()
    # For each file, shift mcs_mask values by cumulative track count
    datasets = []
    track_offset = 0
    for i, f in enumerate(files):
        t1 = time.time()
        print(f"Reading file {i+1}/{len(files)}: {f}")
        ds = xr.open_zarr(f, mask_and_scale=False, decode_times=False)
        # Only shift mcs_mask, not ccs_mask
        mcs_mask = ds["mcs_mask"].copy()
        # Only shift nonzero (background=0)
        mask_nonzero = mcs_mask > 0
        mcs_mask_shifted = mcs_mask.where(~mask_nonzero, mcs_mask + track_offset)
        ds["mcs_mask"] = mcs_mask_shifted
        datasets.append(ds)
        print(f"  track_offset now: {track_offset}")
        # Find max track number in this file for offset
        max_track = int(mcs_mask.max().values)
        track_offset += max_track
        print(f"  max_track in file: {max_track}, new track_offset: {track_offset}")
        print(f"  Done reading file {i+1} ({time.time() - t1:.2f} s elapsed)")

    # Concatenate along time, preserve chunking/compression
    print(f"Concatenating {len(datasets)} datasets...")
    t2 = time.time()
    combined = xr.concat(datasets, dim="time")
    print(f"Concatenation done ({time.time() - t2:.2f} s)")

    # Remove duplicate times, keep last occurrence
    time_vals = combined['time'].values
    _, index = np.unique(time_vals[::-1], return_index=True)
    keep = combined.sizes['time'] - 1 - index
    if len(keep) < combined.sizes['time']:
        print(f"Removing {combined.sizes['time'] - len(keep)} duplicate times (keeping last occurrence)...")
        combined = combined.isel(time=np.sort(keep))
        print(f"Removed duplicate times. New time length: {combined.sizes['time']}")

    # After concatenation, rechunk to match the first file's chunking for each dimension
    first_ds = datasets[0]
    chunk_dict = {dim: chunks[0] for dim, chunks in first_ds.chunks.items()}
    if chunk_dict:
        print(f"Rechunking combined dataset to match original chunk sizes: {chunk_dict}")
        combined = combined.chunk(chunk_dict)

    # Write combined dataset to Zarr using Dask parallelism
    t3 = time.time()
    print(f"Writing output to {out_path} with Dask parallelism ...")
    with ProgressBar():
        combined.to_zarr(out_path, mode="w", compute=True)
    print(f"Combined Zarr written: {out_path}")
    print(f"Writing time elapsed: {time.time() - t3:.2f} s")
    print(f"Total time elapsed: {time.time() - t0:.2f} s")

if __name__ == "__main__":
    main()