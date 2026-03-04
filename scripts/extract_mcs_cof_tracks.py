#!/usr/bin/env python3
"""
Extract MCS Co-occurrence Feature (COF) Flags and Partner Track Numbers

This script reads the COF masks zarr file produced by make_cooccurrence_masks.py,
scans each time step to record which (MCS track ID, COF timestamp) pairs appear in
each overlap category, maps those 6-hourly COF flags back onto the 1-hourly MCS
track statistics time dimension, and additionally stores the actual AR and ETC track
IDs that overlap with each MCS at each time step.

Outputs
-------
When --trackstats is provided the script writes three files:

  1. <source>_mcs_cof_tracks_2d.nc
     NetCDF with 2D arrays shaped (tracks, times) — same layout as base_time in
     the track stats file.

     Binary flag variables (int8):
       cof_flag_mcs_ar      : MCS-AR 2-way overlap present in this hourly step
       cof_flag_mcs_etc     : MCS-ETC 2-way overlap present in this hourly step
       cof_flag_mcs_ar_etc  : MCS-AR-ETC 3-way overlap present in this hourly step
       cof_flag_isolated    : MCS isolated (no overlap) in this hourly step
     Fill value = -1 for time steps with no valid base_time (padding).

     Partner track ID variables (int32):
       ar_tracknum   : AR track ID (COF 1-based convention) overlapping this MCS
                       step; 0 = no overlap; -1 = padding (invalid base_time).
                       When multiple ARs overlap one MCS, the minimum ID is stored.
       etc_tracknum  : ETC track ID overlapping this MCS step; same fill convention.
       ar_tracknum_3way  : AR track ID in the MCS-AR-ETC 3-way overlap only.
       etc_tracknum_3way : ETC track ID in the MCS-AR-ETC 3-way overlap only.

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
COF mask values start at 1 (0 = no MCS / AR / ETC).
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
python extract_mcs_cof_tracks.py --source icon_d3hp003

# Full output with 2D arrays + merged track stats:
python extract_mcs_cof_tracks.py \\
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

# Partner masks: contain AR / ETC track IDs at cells that spatially co-locate
# with an MCS.  Each tuple is (mcs_ref_mask, partner_mask).
# 2-way MCS-AR:        mcs_ar_overlap_mask[cell]  = MCS ID,
#                      ar_mcs_overlap_mask[cell]  = AR ID
# 2-way MCS-ETC:       mcs_etc_overlap_mask[cell] = MCS ID,
#                      etc_mcs_overlap_mask[cell] = ETC ID
# 3-way MCS-AR-ETC:    mcs_ar_etc_overlap_mask[cell]  = MCS ID,
#                      ar_mcs_etc_overlap_mask[cell]  = AR ID,
#                      etc_mcs_ar_overlap_mask[cell]  = ETC ID
_AR_PARTNER_PAIRS = [
    ('mcs_ar_overlap_mask',     'ar_mcs_overlap_mask'),    # 2-way
    ('mcs_ar_etc_overlap_mask', 'ar_mcs_etc_overlap_mask'), # 3-way
]
_ETC_PARTNER_PAIRS = [
    ('mcs_etc_overlap_mask',    'etc_mcs_overlap_mask'),    # 2-way
    ('mcs_ar_etc_overlap_mask', 'etc_mcs_ar_overlap_mask'), # 3-way
]
# 3-way only (for separate _3way output variables)
_AR_PARTNER_PAIRS_3WAY  = [('mcs_ar_etc_overlap_mask', 'ar_mcs_etc_overlap_mask')]
_ETC_PARTNER_PAIRS_3WAY = [('mcs_ar_etc_overlap_mask', 'etc_mcs_ar_overlap_mask')]

# All partner mask variables that need to be loaded from the zarr
_ALL_PARTNER_VARS = sorted({
    v for pairs in (_AR_PARTNER_PAIRS + _ETC_PARTNER_PAIRS)
    for v in pairs
})


# ---------------------------------------------------------------------------
# Step 1 — Scan COF zarr
# ---------------------------------------------------------------------------

def scan_cof_zarr(cof_zarr_path: str, chunk_size: int = 100,
                  logger: logging.Logger = None) -> tuple:
    """
    Scan the COF masks zarr and record, per flag type, the set of COF timestamps
    at which each MCS track ID had an overlap.  Also record the AR / ETC partner
    track IDs associated with each MCS track at each timestamp.

    Parameters
    ----------
    cof_zarr_path : str
        Path to the COF masks zarr store.
    chunk_size : int
        Time steps loaded per chunk.
    logger : logging.Logger, optional

    Returns
    -------
    track_times : dict
        ``{flag_name: {cof_track_id (int): set of int64 timestamps (ns)}}``
    mcs_to_ar : dict
        ``{mcs_cof_id (int): {ts_ns (int): set of ar_cof_id (int)}}``
        Combines 2-way MCS-AR and 3-way MCS-AR-ETC co-occurrences.
    mcs_to_etc : dict
        ``{mcs_cof_id (int): {ts_ns (int): set of etc_cof_id (int)}}``
        Combines 2-way MCS-ETC and 3-way MCS-AR-ETC co-occurrences.
    mcs_to_ar_3way : dict
        Like ``mcs_to_ar`` but restricted to 3-way MCS-AR-ETC cells only.
    mcs_to_etc_3way : dict
        Like ``mcs_to_etc`` but restricted to 3-way MCS-AR-ETC cells only.
    """
    log = logger or logging.getLogger(__name__)
    log.info(f"Opening COF zarr: {cof_zarr_path}")

    ds = xr.open_dataset(cof_zarr_path, engine='zarr', chunks={'time': chunk_size})
    n_times = len(ds.time)
    log.info(f"  {n_times} time steps | cell dim: {ds.sizes.get('cell', '?')}")

    # Check required variables
    required_flag_vars = list(_FLAG_TO_MASK.values())
    missing_flag = [v for v in required_flag_vars if v not in ds.data_vars]
    if missing_flag:
        raise ValueError(f"COF zarr is missing required flag variables: {missing_flag}")

    # Partner variables may not exist in older COF versions — warn rather than crash
    available_partner_vars = [v for v in _ALL_PARTNER_VARS if v in ds.data_vars]
    missing_partner = [v for v in _ALL_PARTNER_VARS if v not in ds.data_vars]
    if missing_partner:
        log.warning(f"Partner mask variables not found (partner track output will be "
                    f"skipped): {missing_partner}")

    # {flag: {track_id: set_of_timestamps_ns}}
    # Store timestamps as int64 nanoseconds for fast set lookup
    track_times: dict[str, dict[int, set]] = {f: {} for f in _FLAG_COLS}

    # {mcs_cof_id: {ts_ns: set_of_partner_ids}}
    mcs_to_ar:       dict[int, dict[int, set]] = {}
    mcs_to_etc:      dict[int, dict[int, set]] = {}
    mcs_to_ar_3way:  dict[int, dict[int, set]] = {}
    mcs_to_etc_3way: dict[int, dict[int, set]] = {}

    vars_to_load = required_flag_vars + available_partner_vars

    n_chunks = int(np.ceil(n_times / chunk_size))
    log.info(f"Scanning {n_chunks} chunk(s) of ≤{chunk_size} time steps...")

    for chunk_idx in range(n_chunks):
        t_start = chunk_idx * chunk_size
        t_end   = min(t_start + chunk_size, n_times)
        log.info(f"  Chunk {chunk_idx + 1}/{n_chunks}: steps {t_start}–{t_end - 1}")

        chunk_ds = ds.isel(time=slice(t_start, t_end))
        chunk_ds = chunk_ds[vars_to_load].compute()

        # COF timestamps for this chunk as int64 nanoseconds
        chunk_times_ns = chunk_ds.time.values.astype('datetime64[ns]').astype(np.int64)

        # ---- Binary flag arrays ----
        for flag, var_name in _FLAG_TO_MASK.items():
            mask_chunk = chunk_ds[var_name].values  # (chunk_t, cell)
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

        # ---- Partner track ID arrays ----
        if available_partner_vars:
            _accumulate_partner(
                chunk_ds, chunk_times_ns,
                _AR_PARTNER_PAIRS, mcs_to_ar,
            )
            _accumulate_partner(
                chunk_ds, chunk_times_ns,
                _ETC_PARTNER_PAIRS, mcs_to_etc,
            )
            _accumulate_partner(
                chunk_ds, chunk_times_ns,
                _AR_PARTNER_PAIRS_3WAY, mcs_to_ar_3way,
            )
            _accumulate_partner(
                chunk_ds, chunk_times_ns,
                _ETC_PARTNER_PAIRS_3WAY, mcs_to_etc_3way,
            )

        del chunk_ds

    ds.close()
    log.info("COF scan complete.")

    for flag in _FLAG_COLS:
        log.info(f"  {flag}: {len(track_times[flag])} unique COF MCS track IDs")
    log.info(f"  AR partner records:  {len(mcs_to_ar)} unique MCS IDs")
    log.info(f"  ETC partner records: {len(mcs_to_etc)} unique MCS IDs")

    return track_times, mcs_to_ar, mcs_to_etc, mcs_to_ar_3way, mcs_to_etc_3way


def _accumulate_partner(chunk_ds, chunk_times_ns: np.ndarray,
                        pairs: list, store: dict):
    """
    For each (mcs_ref_var, partner_var) pair, accumulate partner IDs into store.

    store structure: {mcs_cof_id: {ts_ns: set_of_partner_ids}}
    """
    for mcs_var, partner_var in pairs:
        if mcs_var not in chunk_ds or partner_var not in chunk_ds:
            continue
        mcs_chunk     = chunk_ds[mcs_var].values     # (chunk_t, cell)
        partner_chunk = chunk_ds[partner_var].values  # (chunk_t, cell)

        for ti, ts_ns in enumerate(chunk_times_ns):
            mcs_row     = mcs_chunk[ti].ravel()
            partner_row = partner_chunk[ti].ravel()

            # Only cells where both MCS and partner IDs are positive
            valid = (mcs_row > 0) & (partner_row > 0)
            if not valid.any():
                continue

            mcs_ids     = mcs_row[valid]
            partner_ids = partner_row[valid]

            for mcs_id, partner_id in zip(mcs_ids, partner_ids):
                mcs_int     = int(mcs_id)
                partner_int = int(partner_id)
                if mcs_int not in store:
                    store[mcs_int] = {}
                if ts_ns not in store[mcs_int]:
                    store[mcs_int][ts_ns] = set()
                store[mcs_int][ts_ns].add(partner_int)


# ---------------------------------------------------------------------------
# Step 2 — Build 2D flag arrays aligned to base_time(tracks, times)
# ---------------------------------------------------------------------------

def build_2d_flags(base_time_s: np.ndarray, track_times: dict,
                   cof_window_h: int = 6,
                   logger: logging.Logger = None) -> tuple:
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
        Width of one COF time window in hours (default 6).
    logger : logging.Logger, optional

    Returns
    -------
    flag_2d : dict
        ``{flag: np.ndarray shape (n_tracks, n_times), dtype int8}``
        Values: 1 = overlap present, 0 = no overlap, -1 = padding (invalid base_time).
    valid_mask : np.ndarray, shape (n_tracks, n_times), dtype bool
    base_floored_ns : np.ndarray, shape (n_tracks, n_times), dtype int64
        Each valid cell floored to the nearest COF window boundary (nanoseconds).
        Zero for padding cells.
    """
    log = logger or logging.getLogger(__name__)
    n_tracks, n_times = base_time_s.shape
    window_ns  = int(cof_window_h * 3600 * 1e9)

    log.info(f"Building 2D flag arrays: {n_tracks} tracks × {n_times} times "
             f"(COF window = {cof_window_h} h)")

    # Valid mask: exclude -9999 fill and NaN
    valid_mask = (base_time_s > 0) & ~np.isnan(base_time_s)

    n_valid   = int(valid_mask.sum())
    n_invalid = int((~valid_mask).sum())
    log.info(f"  base_time: {n_valid} valid cells, {n_invalid} fill/NaT cells")

    # Convert valid times to int64 nanoseconds
    base_ns = np.zeros(base_time_s.shape, dtype=np.int64)
    base_ns[valid_mask] = (base_time_s[valid_mask] * 1e9).astype(np.int64)

    # Floor each valid step to its COF window start
    base_floored_ns = (base_ns // window_ns) * window_ns  # (n_tracks, n_times)

    results: dict[str, np.ndarray] = {}

    for flag in _FLAG_COLS:
        tt = track_times[flag]   # {cof_track_id: set_of_ts_ns}
        arr = np.full((n_tracks, n_times), -1, dtype=np.int8)

        for track_idx in range(n_tracks):
            cof_id    = track_idx + 1
            valid_row = valid_mask[track_idx]
            arr[track_idx, valid_row] = 0  # valid but no overlap yet

            if cof_id not in tt:
                continue

            cof_ts_set  = tt[cof_id]
            floored_row = base_floored_ns[track_idx]
            match = np.isin(floored_row[valid_row], list(cof_ts_set))
            arr[track_idx, valid_row] = np.where(match, np.int8(1), np.int8(0))

        results[flag] = arr
        n_overlap = int((arr == 1).sum())
        log.info(f"  {flag}: {n_overlap} / {n_valid} valid step-cells flagged")

    return results, valid_mask, base_floored_ns


# ---------------------------------------------------------------------------
# Step 3 — Build 2D partner track ID arrays
# ---------------------------------------------------------------------------

def build_2d_partner_tracks(base_time_s: np.ndarray,
                             valid_mask: np.ndarray,
                             base_floored_ns: np.ndarray,
                             mcs_to_ar: dict,
                             mcs_to_etc: dict,
                             mcs_to_ar_3way: dict,
                             mcs_to_etc_3way: dict,
                             logger: logging.Logger = None) -> dict:
    """
    Build 2D (tracks × times) arrays containing the AR and ETC track IDs that
    co-occur with each MCS at each time step.

    Parameters
    ----------
    base_time_s : np.ndarray, shape (n_tracks, n_times)
        MCS track stats base_time (epoch seconds).
    valid_mask : np.ndarray, shape (n_tracks, n_times), dtype bool
        True for cells that have a valid base_time (not padding).
    base_floored_ns : np.ndarray, shape (n_tracks, n_times), dtype int64
        Each valid cell floored to the nearest COF window boundary (nanoseconds),
        as returned by ``build_2d_flags``.  Zero for padding cells.
    mcs_to_ar : dict
        ``{mcs_cof_id: {ts_ns: set_of_ar_ids}}`` — combined 2-way + 3-way.
    mcs_to_etc : dict
        ``{mcs_cof_id: {ts_ns: set_of_etc_ids}}`` — combined 2-way + 3-way.
    mcs_to_ar_3way : dict
        ``{mcs_cof_id: {ts_ns: set_of_ar_ids}}`` — 3-way MCS-AR-ETC only.
    mcs_to_etc_3way : dict
        ``{mcs_cof_id: {ts_ns: set_of_etc_ids}}`` — 3-way MCS-AR-ETC only.
    logger : logging.Logger, optional

    Returns
    -------
    dict with keys:
        'ar'        : int32 array (tracks, times) — all AR co-occurrences
        'etc'       : int32 array (tracks, times) — all ETC co-occurrences
        'ar_3way'   : int32 array (tracks, times) — AR in 3-way only
        'etc_3way'  : int32 array (tracks, times) — ETC in 3-way only

    Value convention for all arrays:
        >0    : partner track ID (COF 1-based); minimum ID when multiple partners
        0     : no partner overlap at this valid time step
        -1    : padding (base_time is fill / NaT)
    """
    log = logger or logging.getLogger(__name__)
    n_tracks, n_times = base_time_s.shape

    log.info(f"Building 2D partner track ID arrays: {n_tracks} tracks × {n_times} times")

    configs = [
        ('ar',       mcs_to_ar),
        ('etc',      mcs_to_etc),
        ('ar_3way',  mcs_to_ar_3way),
        ('etc_3way', mcs_to_etc_3way),
    ]

    results: dict[str, np.ndarray] = {}

    for key, store in configs:
        arr = np.full((n_tracks, n_times), -1, dtype=np.int32)

        for track_idx in range(n_tracks):
            cof_id    = track_idx + 1
            valid_row = valid_mask[track_idx]           # (n_times,) bool
            arr[track_idx, valid_row] = 0               # valid, no partner yet

            if cof_id not in store:
                continue

            ts_to_partners = store[cof_id]              # {ts_ns: set_of_ids}
            floored_row    = base_floored_ns[track_idx] # (n_times,)

            # Vectorise: get the set of floored timestamps with any partner
            partner_ts_set = set(ts_to_partners.keys())
            valid_indices  = np.where(valid_row)[0]    # indices of valid times

            for vi in valid_indices:
                ts_ns = int(floored_row[vi])
                if ts_ns in ts_to_partners:
                    # Store minimum partner ID for determinism when multiple overlap
                    arr[track_idx, vi] = min(ts_to_partners[ts_ns])

        n_with_partner = int((arr > 0).sum())
        n_valid_total  = int(valid_mask.sum())
        log.info(f"  {key}: {n_with_partner} / {n_valid_total} valid step-cells "
                 f"have a partner ID")
        results[key] = arr

    return results


# ---------------------------------------------------------------------------
# Step 4 — Compute scalar summaries (any / fraction)
# ---------------------------------------------------------------------------

def compute_scalar_summaries(flag_2d: dict, valid_mask: np.ndarray,
                              logger: logging.Logger = None) -> pd.DataFrame:
    """
    Derive per-track scalar flag and fraction columns from 2D arrays.

    Parameters
    ----------
    flag_2d : dict
        Output of ``build_2d_flags`` (the flag_2d portion).
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
    valid_counts = valid_mask.sum(axis=1).astype(float)  # (n_tracks,)

    rows: dict[str, np.ndarray] = {}
    for flag, frac_col in zip(_FLAG_COLS, _FRAC_COLS):
        arr = flag_2d[flag]           # int8: 1, 0, -1
        overlap = (arr == 1)          # bool
        any_overlap  = overlap.any(axis=1)
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
# Step 5 — Save 2D arrays as netCDF
# ---------------------------------------------------------------------------

def save_2d_netcdf(flag_2d: dict, partner_2d: dict,
                   base_time_s: np.ndarray,
                   n_tracks: int, n_times: int,
                   output_path: str,
                   logger: logging.Logger = None):
    """
    Write the 2D flag and partner track ID arrays to a netCDF file with
    (tracks, times) dimensions that mirrors the MCS track statistics file layout.

    Parameters
    ----------
    flag_2d : dict
        Output of ``build_2d_flags`` (flag arrays, int8).
    partner_2d : dict
        Output of ``build_2d_partner_tracks`` (partner ID arrays, int32).
        May be an empty dict if partner variables were unavailable in the zarr.
    base_time_s : np.ndarray, shape (n_tracks, n_times)
        Original base_time array (epoch seconds) — stored for reference.
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

    # ---- Binary flags ----
    flag_long_names = {
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
                'long_name':     flag_long_names[flag],
                'flag_values':   '1, 0, -1',
                'flag_meanings': '1=overlap 0=no_overlap -1=padding(invalid_time)',
                '_FillValue':    np.int8(-1),
            }
        )

    # ---- Partner track IDs ----
    partner_long_names = {
        'ar':       ('ar_tracknum',      'AR track ID co-occurring with MCS (2-way + 3-way)'),
        'etc':      ('etc_tracknum',     'ETC track ID co-occurring with MCS (2-way + 3-way)'),
        'ar_3way':  ('ar_tracknum_3way', 'AR track ID in MCS-AR-ETC 3-way co-occurrence'),
        'etc_3way': ('etc_tracknum_3way','ETC track ID in MCS-AR-ETC 3-way co-occurrence'),
    }
    if partner_2d:
        for key, (nc_var, long_name) in partner_long_names.items():
            if key not in partner_2d:
                continue
            data_vars[nc_var] = xr.DataArray(
                partner_2d[key],
                dims=['tracks', 'times'],
                attrs={
                    'long_name':   long_name,
                    'flag_values': '>0=partner_track_id, 0=no_overlap, -1=padding',
                    'comment':     ('COF track IDs are 1-based.  When multiple '
                                    'partners overlap one MCS cell, the minimum '
                                    'track ID is stored.'),
                    '_FillValue':  np.int32(-1),
                }
            )

    # ---- base_time for reference ----
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
        'MCS co-occurrence feature (COF) flags and partner AR/ETC track IDs '
        'mapped onto the MCS track statistics (tracks × times) grid.  '
        'Flag values: 1=overlap, 0=no overlap, -1=padding.  '
        'Partner track ID values: >0=track_id, 0=no overlap, -1=padding.'
    )

    # Minimal compression: fast to write, ~3-5× smaller on disk
    encoding = {
        var: {'zlib': True, 'complevel': 1}
        for var in ds_out.data_vars
    }
    ds_out.to_netcdf(output_path, encoding=encoding)
    log.info(f"2D flag + partner tracks netCDF saved → {output_path}")


# ---------------------------------------------------------------------------
# Step 6 — Merge scalar summaries into MCS track statistics DataFrame
# ---------------------------------------------------------------------------

def load_trackstats(trackstats_path: str, logger: logging.Logger = None):
    """
    Load MCS track statistics and return an xarray.Dataset.

    Opens with ``decode_times=False`` and ``mask_and_scale=False`` so that
    ``base_time`` is returned as raw float64 epoch seconds (fill = -9999),
    avoiding xarray silently converting it to datetime64 nanoseconds.
    """
    log = logger or logging.getLogger(__name__)
    ext = Path(trackstats_path).suffix.lower()
    log.info(f"Loading MCS track statistics: {trackstats_path}")

    if ext in ('.nc', '.netcdf', '.nc4'):
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

    Only variables with a single 'tracks' dimension are included.
    """
    log = logger or logging.getLogger(__name__)

    scalar_vars = [v for v in trackstats_ds.data_vars
                   if trackstats_ds[v].dims == ('tracks',)]
    log.info(f"  1D track variables to include: {scalar_vars}")

    stats_df = trackstats_ds[scalar_vars].to_dataframe().reset_index()

    summary = scalar_df.reset_index()  # 'tracks' becomes a column
    merged  = stats_df.merge(summary, on='tracks', how='left')

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
        description="Extract MCS COF flags and AR/ETC partner track IDs from COF masks zarr."
    )
    parser.add_argument('--source', required=True,
                        help='Source name, used to locate '
                             '<cof-dir>/<source>_cofmasks_hp8_v1.zarr '
                             '(e.g. icon_d3hp003, scream, IMERGv7)')
    parser.add_argument('--trackstats', default=None,
                        help='Path to MCS track statistics netCDF (.nc).  '
                             'Required for 2D output and fraction columns.')
    parser.add_argument('--cof-dir', default=None,
                        help='Directory containing COF masks zarr stores '
                             '(default: /pscratch/sd/w/wcmca1/hackathon/cof_masks)')
    parser.add_argument('--output-dir', default=None,
                        help='Output directory '
                             '(default: <cof-dir>/stats)')
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
    tracks_2d_nc   = os.path.join(out_dir, f"{source_name}_mcs_cof_tracks_2d.nc")
    merged_parquet = os.path.join(out_dir, f"{source_name}_mcs_trackstats_cof.parquet")

    logger.info("=" * 70)
    logger.info(f"MCS COF FLAG + PARTNER TRACK EXTRACTION  |  source: {source_name}")
    logger.info("=" * 70)

    # ------------------------------------------------------------------
    # Step 1: Scan COF zarr → per-track-per-time overlap records
    # ------------------------------------------------------------------
    logger.info("Step 1: Scanning COF zarr for track-time overlap records...")
    track_times, mcs_to_ar, mcs_to_etc, mcs_to_ar_3way, mcs_to_etc_3way = scan_cof_zarr(
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

        base_time_s = ts_ds['base_time'].values.astype(float)
        n_tracks, n_times = base_time_s.shape

        logger.info("Step 3: Building 2D flag arrays (tracks × times)...")
        flag_2d, valid_mask, base_floored_ns = build_2d_flags(
            base_time_s=base_time_s,
            track_times=track_times,
            cof_window_h=args.cof_window,
            logger=logger,
        )

        logger.info("Step 4: Building 2D partner track ID arrays (tracks × times)...")
        partner_2d = build_2d_partner_tracks(
            base_time_s=base_time_s,
            valid_mask=valid_mask,
            base_floored_ns=base_floored_ns,
            mcs_to_ar=mcs_to_ar,
            mcs_to_etc=mcs_to_etc,
            mcs_to_ar_3way=mcs_to_ar_3way,
            mcs_to_etc_3way=mcs_to_etc_3way,
            logger=logger,
        )

        logger.info("Step 5: Saving 2D arrays to netCDF...")
        save_2d_netcdf(
            flag_2d=flag_2d,
            partner_2d=partner_2d,
            base_time_s=base_time_s,
            n_tracks=n_tracks,
            n_times=n_times,
            output_path=tracks_2d_nc,
            logger=logger,
        )

        logger.info("Step 6: Computing scalar summaries (any / fraction)...")
        scalar_df = compute_scalar_summaries(
            flag_2d=flag_2d,
            valid_mask=valid_mask,
            logger=logger,
        )

        scalar_df.to_parquet(flags_parquet, index=True, engine='pyarrow')
        logger.info(f"Scalar flags parquet saved → {flags_parquet}")

        logger.info("Step 7: Merging scalars into track statistics parquet...")
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
        logger.info(f"  2D tracks netCDF   : {tracks_2d_nc}")
        logger.info(f"  Merged track stats : {merged_parquet}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
