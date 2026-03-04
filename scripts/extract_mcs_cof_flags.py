#!/usr/bin/env python3
"""
Extract MCS Co-occurrence Feature (COF) Flags

This script reads the COF masks zarr file produced by make_cooccurrence_masks.py,
scans each time step to record which (MCS track ID, COF timestamp) pairs appear in
each overlap category, then maps those 6-hourly COF flags back onto the 1-hourly
MCS track statistics time dimension.

Outputs
-------
When --trackstats is provided the script writes three files:

  1. <source>_mcs_cof_flags_2d.nc
     NetCDF with 2D flag arrays shaped (tracks, times) — same layout as base_time in
     the track stats file.  A cell is 1 if the corresponding 1-hourly step falls inside
     a 6-hourly COF window that had an overlap for that track.
       cof_flag_mcs_ar      : MCS-AR 2-way overlap present in this hourly step
       cof_flag_mcs_etc     : MCS-ETC 2-way overlap present in this hourly step
       cof_flag_mcs_ar_etc  : MCS-AR-ETC 3-way overlap present in this hourly step
       cof_flag_isolated    : MCS isolated (no overlap) in this hourly step
     Fill value = -1 for time steps that have no valid base_time (padding).

  2. <source>_mcs_cof_flags.parquet
     One row per MCS track (indexed 0-based, matching the tracks dimension).
     Scalar summary columns:
       flag_mcs_ar / flag_mcs_etc / flag_mcs_ar_etc / flag_isolated
           True if the track has ≥1 overlapping hourly step (any time)
       frac_mcs_ar / frac_mcs_etc / frac_mcs_ar_etc / frac_isolated
           Fraction of the track's valid lifetime steps that have overlap

  3. <source>_mcs_trackstats_cof.parquet
     The original MCS track statistics (from --trackstats), converted to a flat
     DataFrame with the scalar flag/fraction columns appended.

If --trackstats is NOT provided, only a simple per-track flags parquet is written
(Step 1 from the old behaviour), without 2D arrays or fractions.

Track ID convention
-------------------
COF mask values start at 1 (0 = no MCS).
MCS track stats use a 0-based 'tracks' dimension index.
Mapping: stats index i  ↔  COF mask value i+1.

COF time → 1-hourly mapping
----------------------------
make_mcs_swath_masks.py floors input times to the nearest N-hourly boundary
(default N=6) using dt.floor(f'{N}h').  A COF timestamp T therefore covers
all 1-hourly steps in [T, T+N*h).  The window size is set with --cof-window
(default 6).

Usage
-----
# Flags parquet only (no trackstats):
python extract_mcs_cof_flags.py --source icon_d3hp003

# Full output with 2D arrays + merged track stats:
python extract_mcs_cof_flags.py \\
    --source icon_d3hp003 \\
    --trackstats /path/to/mcs_tracks_final.nc \\
    --cof-window 6 \\
    --chunk-size 100

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import os
import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Mapping: scalar flag name → COF mask variable containing MCS track IDs
_FLAG_TO_MASK = {
    'flag_mcs_ar':     'mcs_ar_overlap_mask',
    'flag_mcs_etc':    'mcs_etc_overlap_mask',
    'flag_mcs_ar_etc': 'mcs_ar_etc_overlap_mask',
    'flag_isolated':   'mcs_isolated_mask',
}

_FLAG_COLS   = list(_FLAG_TO_MASK.keys())
_FRAC_COLS   = [f.replace('flag_', 'frac_') for f in _FLAG_COLS]
_NC_VAR_COLS = [f'cof_{f}' for f in _FLAG_COLS]   # names inside the netCDF


# ---------------------------------------------------------------------------
# Step 1 — Scan COF zarr: build {flag: {track_id: set_of_cof_timestamps}}
# ---------------------------------------------------------------------------

def scan_cof_zarr(cof_zarr_path: str, chunk_size: int = 100,
                  logger: logging.Logger = None) -> dict:
    """
    Scan the COF masks zarr and record, per flag type, the set of COF timestamps
    at which each MCS track ID had an overlap.

    Parameters
    ----------
    cof_zarr_path : str
        Path to the COF masks zarr store.
    chunk_size : int
        Time steps loaded per chunk.
    logger : logging.Logger, optional

    Returns
    -------
    dict
        ``{flag_name: {cof_track_id (int): set of np.datetime64 timestamps}}``
        ``cof_track_id`` uses the COF convention (starts at 1).
    """
    log = logger or logging.getLogger(__name__)
    log.info(f"Opening COF zarr: {cof_zarr_path}")

    ds = xr.open_dataset(cof_zarr_path, engine='zarr', chunks={'time': chunk_size})
    n_times = len(ds.time)
    log.info(f"  {n_times} time steps | cell dim: {ds.sizes.get('cell', '?')}")

    missing = [v for v in _FLAG_TO_MASK.values() if v not in ds.data_vars]
    if missing:
        raise ValueError(f"COF zarr is missing required variables: {missing}")

    # {flag: {track_id: set_of_timestamps_ns}}
    # Store timestamps as int64 nanoseconds for fast set lookup
    track_times: dict[str, dict[int, set]] = {f: {} for f in _FLAG_COLS}

    n_chunks = int(np.ceil(n_times / chunk_size))
    log.info(f"Scanning {n_chunks} chunk(s) of ≤{chunk_size} time steps...")

    for chunk_idx in range(n_chunks):
        t_start = chunk_idx * chunk_size
        t_end   = min(t_start + chunk_size, n_times)
        log.info(f"  Chunk {chunk_idx + 1}/{n_chunks}: steps {t_start}–{t_end - 1}")

        # Load only the mask variables and their times
        chunk_ds = ds.isel(time=slice(t_start, t_end))
        chunk_ds = chunk_ds[list(_FLAG_TO_MASK.values())].compute()

        # COF timestamps for this chunk as int64 nanoseconds
        chunk_times_ns = chunk_ds.time.values.astype('datetime64[ns]').astype(np.int64)

        for flag, var_name in _FLAG_TO_MASK.items():
            mask_chunk = chunk_ds[var_name].values  # shape: (chunk_t, cell)
            tt = track_times[flag]

            for ti, ts_ns in enumerate(chunk_times_ns):
                row = mask_chunk[ti].ravel()
                ids = row[row > 0]
                if ids.size == 0:
                    continue
                for tid in np.unique(ids):
                    tid_int = int(tid)
                    if tid_int not in tt:
                        tt[tid_int] = set()
                    tt[tid_int].add(ts_ns)

        del chunk_ds

    ds.close()
    log.info("COF scan complete.")

    # Summary
    for flag in _FLAG_COLS:
        log.info(f"  {flag}: {len(track_times[flag])} unique COF track IDs")

    return track_times


# ---------------------------------------------------------------------------
# Step 2 — Build 2D flag arrays aligned to base_time(tracks, times)
# ---------------------------------------------------------------------------

def build_2d_flags(base_time_s: np.ndarray, track_times: dict,
                   cof_window_h: int = 6,
                   logger: logging.Logger = None) -> dict:
    """
    Map the 6-hourly COF overlap records onto the 2D (tracks, times) grid of
    the MCS track statistics file.

    Parameters
    ----------
    base_time_s : np.ndarray, shape (n_tracks, n_times), dtype float64
        Epoch seconds from the MCS track stats ``base_time`` variable.
        Fill value is -9999 (invalid / padding steps).
    track_times : dict
        Output of ``scan_cof_zarr``:
        ``{flag: {cof_track_id: set_of_timestamps_ns (int64)}}``.
        COF track IDs start at 1; stats track indices are 0-based.
    cof_window_h : int
        Width of one COF time window in hours (default 6, matching swath
        aggregation in make_mcs_swath_masks.py).
    logger : logging.Logger, optional

    Returns
    -------
    dict
        ``{flag: np.ndarray shape (n_tracks, n_times), dtype int8}``
        Values: 1 = overlap present, 0 = no overlap, -1 = padding (invalid base_time).
    """
    log = logger or logging.getLogger(__name__)
    n_tracks, n_times = base_time_s.shape
    window_ns  = int(cof_window_h * 3600 * 1e9)  # window width in nanoseconds

    log.info(f"Building 2D flag arrays: {n_tracks} tracks × {n_times} times "
             f"(COF window = {cof_window_h} h)")

    # Build valid mask: exclude fill values (-9999) AND any NaN that may
    # arise from upstream processing.  Valid epoch seconds for 2019-2021
    # are ~1.5e9, so base_time_s > 0 is a robust positive-time guard that
    # also rejects -9999 fill in one comparison.
    valid_mask = (base_time_s > 0) & ~np.isnan(base_time_s)

    n_valid   = int(valid_mask.sum())
    n_invalid = int((~valid_mask).sum())
    log.info(f"  base_time: {n_valid} valid cells, {n_invalid} fill/NaT cells")

    # Convert valid times to int64 nanoseconds; only cast valid elements to
    # avoid the NaN/fill→int64 overflow warning.
    base_ns = np.zeros(base_time_s.shape, dtype=np.int64)
    base_ns[valid_mask] = (base_time_s[valid_mask] * 1e9).astype(np.int64)

    # For each valid time step, the corresponding COF window start is:
    #   floor(base_ns / window_ns) * window_ns
    base_floored_ns = (base_ns // window_ns) * window_ns  # shape (n_tracks, n_times)

    results: dict[str, np.ndarray] = {}

    for flag in _FLAG_COLS:
        tt = track_times[flag]   # {cof_track_id: set_of_ts_ns}
        arr = np.full((n_tracks, n_times), -1, dtype=np.int8)  # default = padding

        for track_idx in range(n_tracks):
            # COF mask IDs are 1-based; stats dimension is 0-based
            cof_id = track_idx + 1
            valid_row = valid_mask[track_idx]            # shape (n_times,)
            arr[track_idx, valid_row] = 0                # valid but no overlap yet

            if cof_id not in tt:
                continue  # this track never appeared in this COF category

            cof_ts_set = tt[cof_id]                      # set of int64 ns values

            floored_row = base_floored_ns[track_idx]     # shape (n_times,)
            # Vectorized check: which floored times are in the COF set?
            match = np.isin(floored_row[valid_row], list(cof_ts_set))
            arr[track_idx, valid_row] = np.where(match, np.int8(1), np.int8(0))

        results[flag] = arr
        n_overlap = int((arr == 1).sum())
        n_valid   = int(valid_mask.sum())
        log.info(f"  {flag}: {n_overlap} / {n_valid} valid step-cells flagged")

    return results


# ---------------------------------------------------------------------------
# Step 3 — Compute scalar summaries (any / fraction)
# ---------------------------------------------------------------------------

def compute_scalar_summaries(flag_2d: dict, valid_mask: np.ndarray,
                              logger: logging.Logger = None) -> pd.DataFrame:
    """
    Derive per-track scalar flag and fraction columns from 2D arrays.

    Parameters
    ----------
    flag_2d : dict
        Output of ``build_2d_flags``.
    valid_mask : np.ndarray, shape (n_tracks, n_times), dtype bool
        True where base_time is valid (not padding).
    logger : logging.Logger, optional

    Returns
    -------
    pd.DataFrame
        Index = track index (0-based).  Columns:
        ``flag_mcs_ar``, ``flag_mcs_etc``, ``flag_mcs_ar_etc``, ``flag_isolated``,
        ``frac_mcs_ar``, ``frac_mcs_etc``, ``frac_mcs_ar_etc``, ``frac_isolated``.
    """
    log = logger or logging.getLogger(__name__)
    n_tracks = next(iter(flag_2d.values())).shape[0]
    valid_counts = valid_mask.sum(axis=1).astype(float)  # (n_tracks,)

    rows: dict[str, np.ndarray] = {}
    for flag, frac_col in zip(_FLAG_COLS, _FRAC_COLS):
        arr = flag_2d[flag]           # int8: 1, 0, -1
        overlap = (arr == 1)          # bool
        any_overlap  = overlap.any(axis=1)                                  # bool (n_tracks,)
        frac_overlap = np.where(valid_counts > 0,
                                overlap.sum(axis=1) / valid_counts,
                                0.0)
        rows[flag]     = any_overlap
        rows[frac_col] = frac_overlap

    df = pd.DataFrame(rows)
    df.index.name = 'tracks'

    log.info("Scalar summary (# tracks with any overlap):")
    for flag in _FLAG_COLS:
        log.info(f"  {flag}: {df[flag].sum()} "
                 f"| mean frac: {df[flag.replace('flag_','frac_')].mean():.3f}")

    return df


# ---------------------------------------------------------------------------
# Step 4 — Save 2D flags as netCDF
# ---------------------------------------------------------------------------

def save_2d_flags_netcdf(flag_2d: dict, base_time_s: np.ndarray,
                         n_tracks: int, n_times: int,
                         output_path: str,
                         logger: logging.Logger = None):
    """
    Write the 2D flag arrays to a netCDF file with (tracks, times) dimensions
    that mirrors the MCS track statistics file layout.

    Parameters
    ----------
    flag_2d : dict
        Output of ``build_2d_flags``.
    base_time_s : np.ndarray, shape (n_tracks, n_times)
        Original base_time array (epoch seconds) — stored as a coordinate.
    n_tracks, n_times : int
        Dimension sizes.
    output_path : str
        Destination file path.
    logger : logging.Logger, optional
    """
    log = logger or logging.getLogger(__name__)

    coords = {
        'tracks': np.arange(n_tracks, dtype=np.int64),
        'times':  np.arange(n_times,  dtype=np.int64),
    }

    data_vars = {}
    long_names = {
        'flag_mcs_ar':     'MCS-AR 2-way co-occurrence flag',
        'flag_mcs_etc':    'MCS-ETC 2-way co-occurrence flag',
        'flag_mcs_ar_etc': 'MCS-AR-ETC 3-way co-occurrence flag',
        'flag_isolated':   'Isolated MCS flag (no co-occurrence)',
    }
    for flag, nc_var in zip(_FLAG_COLS, _NC_VAR_COLS):
        data_vars[nc_var] = xr.DataArray(
            flag_2d[flag],
            dims=['tracks', 'times'],
            attrs={
                'long_name':   long_names[flag],
                'flag_values': '1, 0, -1',
                'flag_meanings': '1=overlap 0=no_overlap -1=padding(invalid_time)',
                '_FillValue':  np.int8(-1),
            }
        )

    # Also store base_time for reference
    data_vars['base_time'] = xr.DataArray(
        base_time_s,
        dims=['tracks', 'times'],
        attrs={
            'long_name': 'Base time of MCS track step',
            'units': 'Seconds since 1970-1-1',
            '_FillValue': -9999.0,
        }
    )

    ds_out = xr.Dataset(data_vars, coords=coords)
    ds_out.attrs['description'] = (
        'MCS co-occurrence feature (COF) flags mapped onto the MCS track '
        'statistics (tracks × times) grid.  Values: 1=overlap, 0=no overlap, '
        '-1=padding (base_time fill).'
    )

    # Minimal compression (complevel=1): fast to write, ~3-5× smaller on disk
    encoding = {
        var: {'zlib': True, 'complevel': 1}
        for var in ds_out.data_vars
    }
    ds_out.to_netcdf(output_path, encoding=encoding)
    log.info(f"2D flag netCDF saved → {output_path}")


# ---------------------------------------------------------------------------
# Step 5 — Merge scalar summaries into MCS track statistics DataFrame
# ---------------------------------------------------------------------------

def load_trackstats(trackstats_path: str, logger: logging.Logger = None):
    """
    Load MCS track statistics and return (xarray.Dataset, flat pd.DataFrame).

    Returns the dataset so callers can access base_time as a raw array without
    converting the entire file to a DataFrame first.

    Opens with ``decode_times=False`` and ``mask_and_scale=False`` so that
    ``base_time`` is returned as raw float64 epoch seconds (fill = -9999),
    avoiding xarray silently converting it to datetime64 nanoseconds.
    """
    log = logger or logging.getLogger(__name__)
    ext = Path(trackstats_path).suffix.lower()
    log.info(f"Loading MCS track statistics: {trackstats_path}")

    if ext in ('.nc', '.netcdf', '.nc4'):
        # decode_times=False keeps base_time as raw float epoch seconds
        # mask_and_scale=False keeps the -9999 fill value as-is (no NaN)
        ds = xr.open_dataset(trackstats_path, decode_times=False, mask_and_scale=False)
    elif ext == '.parquet':
        raise ValueError("--trackstats must be a netCDF file (.nc) so that "
                         "the 2D base_time array is accessible.")
    else:
        raise ValueError(f"Unsupported format: {ext}.  Expected .nc / .nc4")

    log.info(f"  tracks={ds.sizes['tracks']}, times={ds.sizes['times']}")
    return ds


def merge_scalars_into_trackstats(trackstats_ds: xr.Dataset,
                                   scalar_df: pd.DataFrame,
                                   output_path: str,
                                   logger: logging.Logger = None):
    """
    Flatten the MCS track statistics dataset to a DataFrame, join scalar COF
    summaries, and save as parquet.

    Only variables with a single 'tracks' dimension (i.e., per-track scalars)
    are included to keep the output manageable.  Variables with additional
    dimensions (e.g., base_time[tracks, times]) are excluded.
    """
    log = logger or logging.getLogger(__name__)

    # Collect 1D (tracks-only) variables
    scalar_vars = [v for v in trackstats_ds.data_vars
                   if trackstats_ds[v].dims == ('tracks',)]
    log.info(f"  1D track variables to include: {scalar_vars}")

    stats_df = trackstats_ds[scalar_vars].to_dataframe().reset_index()

    # scalar_df index is 0-based track index; rename to join on 'tracks'
    summary = scalar_df.reset_index()  # 'tracks' becomes a column

    merged = stats_df.merge(summary, on='tracks', how='left')

    # Fill unmatched (shouldn't happen if stats and COF cover same period)
    for col in _FLAG_COLS:
        if col in merged.columns:
            merged[col] = merged[col].fillna(False).astype(bool)
    for col in _FRAC_COLS:
        if col in merged.columns:
            merged[col] = merged[col].fillna(0.0).astype(float)

    merged.to_parquet(output_path, index=False, engine='pyarrow')
    log.info(f"Merged track stats saved → {output_path}")
    log.info(f"  Shape: {merged.shape}")
    return merged


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract MCS co-occurrence feature (COF) flags from COF masks zarr."
    )
    parser.add_argument('--source', required=True,
                        help='Source name, used to locate '
                             '<cof-dir>/<source>_cofmasks_hp8_v1.zarr '
                             '(e.g. icon_d3hp003, scream, IMERGv7)')
    parser.add_argument('--trackstats', default=None,
                        help='Path to MCS track statistics netCDF (.nc).  '
                             'Required for 2D flag output and fraction columns.')
    parser.add_argument('--cof-dir', default=None,
                        help='Directory containing COF masks zarr stores '
                             '(default: /pscratch/sd/w/wcmca1/hackathon/cof_masks)')
    parser.add_argument('--output-dir', default=None,
                        help='Output directory (default: same as --cof-dir)')
    parser.add_argument('--cof-window', type=int, default=6,
                        help='COF aggregation window in hours — must match the '
                             '--aggregation-window used in make_mcs_swath_masks.py '
                             '(default: 6)')
    parser.add_argument('--chunk-size', type=int, default=100,
                        help='Time steps per processing chunk (default: 100)')
    return parser.parse_args()


def main():
    setup_logging()
    logger = logging.getLogger(__name__)
    args = parse_args()

    source_name = args.source
    cof_root    = args.cof_dir or "/pscratch/sd/w/wcmca1/hackathon/cof_masks"
    cof_zarr    = os.path.join(cof_root, f"{source_name}_cofmasks_hp8_v1.zarr")
    out_dir     = args.output_dir or os.path.join(cof_root, "stats")
    os.makedirs(out_dir, exist_ok=True)

    flags_parquet  = os.path.join(out_dir, f"{source_name}_mcs_cof_flags.parquet")
    flags_2d_nc    = os.path.join(out_dir, f"{source_name}_mcs_cof_flags_2d.nc")
    merged_parquet = os.path.join(out_dir, f"{source_name}_mcs_trackstats_cof.parquet")

    logger.info("=" * 70)
    logger.info(f"MCS COF FLAG EXTRACTION  |  source: {source_name}")
    logger.info("=" * 70)

    # ------------------------------------------------------------------
    # Step 1: Scan COF zarr → per-track-per-time overlap records
    # ------------------------------------------------------------------
    logger.info("Step 1: Scanning COF zarr for track-time overlap records...")
    track_times = scan_cof_zarr(
        cof_zarr_path=cof_zarr,
        chunk_size=args.chunk_size,
        logger=logger,
    )

    # ------------------------------------------------------------------
    # Step 2: If trackstats provided, build 2D arrays and save netCDF
    # ------------------------------------------------------------------
    if args.trackstats:
        logger.info("Step 2: Loading MCS track statistics...")
        ts_ds = load_trackstats(args.trackstats, logger=logger)

        # base_time: shape (n_tracks, n_times), epoch seconds, fill=-9999
        base_time_s = ts_ds['base_time'].values.astype(float)
        n_tracks, n_times = base_time_s.shape
        valid_mask = (base_time_s > 0) & ~np.isnan(base_time_s)

        logger.info("Step 3: Building 2D flag arrays (tracks × times)...")
        flag_2d = build_2d_flags(
            base_time_s=base_time_s,
            track_times=track_times,
            cof_window_h=args.cof_window,
            logger=logger,
        )

        logger.info("Step 4: Saving 2D flags to netCDF...")
        save_2d_flags_netcdf(
            flag_2d=flag_2d,
            base_time_s=base_time_s,
            n_tracks=n_tracks,
            n_times=n_times,
            output_path=flags_2d_nc,
            logger=logger,
        )

        logger.info("Step 5: Computing scalar summaries (any / fraction)...")
        scalar_df = compute_scalar_summaries(
            flag_2d=flag_2d,
            valid_mask=valid_mask,
            logger=logger,
        )

        scalar_df.to_parquet(flags_parquet, index=True, engine='pyarrow')
        logger.info(f"Scalar flags parquet saved → {flags_parquet}")

        logger.info("Step 6: Merging scalars into track statistics parquet...")
        merge_scalars_into_trackstats(
            trackstats_ds=ts_ds,
            scalar_df=scalar_df,
            output_path=merged_parquet,
            logger=logger,
        )
        ts_ds.close()

    else:
        # ------------------------------------------------------------------
        # Fallback: no trackstats — write a simple per-track flags parquet
        # (any-overlap only, no fractions or 2D arrays)
        # ------------------------------------------------------------------
        logger.info("Step 2: Building simple per-track flags (no trackstats provided)...")
        all_ids = sorted(
            set(track_times['flag_isolated'])
            | set(track_times['flag_mcs_ar'])
            | set(track_times['flag_mcs_etc'])
            | set(track_times['flag_mcs_ar_etc'])
        )
        if not all_ids:
            logger.warning("No MCS track IDs found in any COF mask variable!")
        else:
            df = pd.DataFrame({'mcs_cof_track_id': np.array(all_ids, dtype=np.int64)})
            for flag in _FLAG_COLS:
                df[flag] = df['mcs_cof_track_id'].isin(track_times[flag])
            df.to_parquet(flags_parquet, index=False, engine='pyarrow')
            logger.info(f"Simple flags parquet saved → {flags_parquet}")
            for flag in _FLAG_COLS:
                logger.info(f"  {flag}: {df[flag].sum()} tracks")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    logger.info("=" * 70)
    logger.info("DONE")
    logger.info(f"  Flags parquet      : {flags_parquet}")
    if args.trackstats:
        logger.info(f"  2D flags netCDF    : {flags_2d_nc}")
        logger.info(f"  Merged track stats : {merged_parquet}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
