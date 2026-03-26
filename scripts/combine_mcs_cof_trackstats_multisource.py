#!/usr/bin/env python3
"""
Combine MCS Track Statistics + COF 2D Overlap Data (Multi-Source)

For each source (OBS/IMERGv7, SCREAM, ICON, UM, NICAM, CASESM2), loads:
  - MCS track statistics netCDF  (from /pscratch/.../mcs/<source>/stats/)
  - COF 2D tracks netCDF         (from /pscratch/.../cof_masks/stats/)

Processing steps
----------------
  1. Load & merge: open track-stats and COF 2D files per source, merge on
     the shared (tracks, times) grid.
  1b. Tropical filter: optionally remove tracks whose lifetime-median
     |meanlat| falls within ±tropics_lat_threshold (default 20°).  The
     filter returns integer keep-indices; the original unmodified dataset
     is passed to the conversion step to avoid slow chunk-by-chunk reads
     on compressed netCDF variables.
  2. Convert: tidy DataFrame (one row per valid MCS time step) with
     dimension reduction:
       (tracks, times, mergers):
         - merge_cloudnumber, split_cloudnumber → discarded
         - merge_ccs_area, split_ccs_area       → summed over mergers dim
       (tracks, times, nmaxpf):
         - all PF variables                     → nmaxpf index 0 (largest PF)
       (tracks, times), (tracks,):              kept as-is
  3. Save: concatenate all sources and write as a snappy-compressed parquet.

Usage
-----
python combine_mcs_cof_trackstats_multisource.py
python combine_mcs_cof_trackstats_multisource.py --sources scream IMERGv7
python combine_mcs_cof_trackstats_multisource.py --tropics-lat-threshold 0  # keep all latitudes
python combine_mcs_cof_trackstats_multisource.py --output-dir /my/output --sources scream IMERGv7

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import os
import time
import logging
import argparse
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
# Configuration
# ---------------------------------------------------------------------------

MCS_DIR     = "/pscratch/sd/w/wcmca1/hackathon/mcs/"
COF_DIR     = "/pscratch/sd/w/wcmca1/hackathon/cof_masks/stats/"
OUTPUT_DIR  = "/pscratch/sd/w/wcmca1/hackathon/cof_masks/stats/"

DEFAULT_SOURCES = [
    "IMERGv7",
    "scream",
    "icon_d3hp003",
    "um_glm_n2560_RAL3p3",
    "nicam_gl11",
    "casesm2_10km_nocumulus",
]

TRACKSTATS_FILES = {
    "IMERGv7":                "mcs_tracks_final_20190101.0000_20220101.0100.nc",
    "scream":                 "mcs_tracks_final_20190801.0000_20200901.0000.nc",
    "icon_d3hp003":           "mcs_tracks_final_20200102.0000_20201231.2330.nc",
    "um_glm_n2560_RAL3p3":   "mcs_tracks_final_20200201.0000_20210301.0000.nc",
    "nicam_gl11":             "mcs_tracks_final_20200301.0000_20210301.0000.nc",
    "casesm2_10km_nocumulus": "mcs_tracks_final_20200301.0000_20210301.0000.nc",
}

SOURCE_DISPLAY_NAMES = {
    "IMERGv7":                "OBS",
    "scream":                 "SCREAM",
    "icon_d3hp003":           "ICON",
    "um_glm_n2560_RAL3p3":   "UM",
    "nicam_gl11":             "NICAM",
    "casesm2_10km_nocumulus": "CASESM2",
}


# ---------------------------------------------------------------------------
# Step 1: Load and merge track stats + COF 2D file per source
# ---------------------------------------------------------------------------

def load_combined(source, mcs_root, cof_stats_dir, trackstats_map,
                  logger=None):
    """
    Open the MCS track statistics netCDF and the COF 2D tracks netCDF for
    one source and return a merged xarray.Dataset.

    Both files are opened with decode_times=False, mask_and_scale=False so
    that base_time stays as raw float64 epoch seconds.
    """
    log = logger or logging.getLogger(__name__)

    ts_path  = os.path.join(mcs_root, source, "stats", trackstats_map[source])
    cof_path = os.path.join(cof_stats_dir, f"{source}_mcs_cof_tracks_2d.nc")

    if not os.path.exists(ts_path):
        raise FileNotFoundError(f"Track stats not found: {ts_path}")
    if not os.path.exists(cof_path):
        raise FileNotFoundError(f"COF 2D file not found: {cof_path}")

    ds_ts  = xr.open_dataset(ts_path,  decode_times=False, mask_and_scale=False)
    ds_cof = xr.open_dataset(cof_path, decode_times=False, mask_and_scale=False)

    # Drop base_time from the COF file — already present in track stats
    cof_vars = [v for v in ds_cof.data_vars if v != 'base_time']
    ds_merged = xr.merge([ds_ts, ds_cof[cof_vars]], join='exact')

    log.info(f"  {source}: {ds_merged.sizes['tracks']} tracks, "
             f"{ds_merged.sizes['times']} times | "
             f"COF vars: {cof_vars}")
    return ds_merged


# ---------------------------------------------------------------------------
# Step 1b: Filter out tropical MCS tracks
# ---------------------------------------------------------------------------

def filter_extratropical_tracks(ds, tropics_lat_threshold=20.0, logger=None):
    """
    Remove MCS tracks whose lifetime-median meanlat falls within the tropics.

    A track is considered tropical if:
        abs(median(meanlat along times)) < tropics_lat_threshold

    Only valid (non-fill) meanlat values are used for the median.
    Fill values (-9999) are masked before computing the median.

    Parameters
    ----------
    ds : xarray.Dataset
        Merged MCS track stats dataset.  Must contain 'meanlat' (tracks, times).
    tropics_lat_threshold : float
        Absolute latitude threshold in degrees (default 20.0).
        Tracks whose median |meanlat| < threshold are removed.
    logger : logging.Logger, optional

    Returns
    -------
    xarray.Dataset
        Dataset with tropical tracks removed.
    """
    log = logger or logging.getLogger(__name__)
    n_before = ds.sizes['tracks']

    meanlat = ds['meanlat'].values.astype(float)   # (tracks, times)
    # Mask fill values (-9999) and any negative-fill sentinel
    meanlat_masked = np.where(meanlat < -900, np.nan, meanlat)

    # Lifetime-median ignoring NaN padding
    median_lat = np.nanmedian(meanlat_masked, axis=1)  # (tracks,)

    # Keep tracks outside the tropics
    keep = np.abs(median_lat) >= tropics_lat_threshold
    keep_indices = np.where(keep)[0]   # integer indices into original ds

    n_after   = keep_indices.size
    n_removed = n_before - n_after
    log.info(f"  Tropical filter (|meanlat_median| < {tropics_lat_threshold}°): "
             f"{n_removed} tracks removed, {n_after} tracks retained "
             f"({n_removed / n_before:.1%} removed)")
    # Return just the keep_indices — caller uses the original dataset for I/O
    # to avoid slow chunk-by-chunk reads on compressed netCDF variables.
    return keep_indices


# ---------------------------------------------------------------------------
# Step 2: Convert merged dataset to tidy Pandas DataFrame
# ---------------------------------------------------------------------------

def ds_to_dataframe(ds, keep_indices=None, time_res_h=1.0, name='', logger=None):
    """
    Convert a merged MCS track stats + COF xarray.Dataset to a tidy DataFrame.

    Only valid time steps (determined by track_duration) are kept — padded
    cells are never materialised in memory.

    Parameters
    ----------
    ds : xarray.Dataset
        The ORIGINAL (pre-filter) merged dataset.  Must not have lazy isel
        fancy-indexing applied — this ensures each ds[v].values call makes
        a single contiguous netCDF4 read (avoids chunk-by-chunk reads on
        compressed variables).
    keep_indices : np.ndarray of int or None
        Integer indices of tracks to retain (output of
        filter_extratropical_tracks).  If None, all tracks are used.

    Dimension handling
    ------------------
    (tracks,)               → broadcast to all valid steps
    (tracks, times)         → fancy-index [orig_track_idx, time_idx]
    (tracks, times, mergers)
        merge_cloudnumber, split_cloudnumber   → discarded (never read)
        merge_ccs_area, split_ccs_area         → summed over mergers → 2D
        all other mergers vars                 → skipped
    (tracks, times, nmaxpf) → slice [:, :, 0] (largest PF) → 2D
    anything else           → skipped

    Extra columns added
    -------------------
    relative_step    : integer step index from track start (0-based)
    relative_time_h  : relative time in hours from track initiation
    track_duration_h : total track lifetime in hours
    dataset          : source label string
    """
    log = logger or logging.getLogger(__name__)

    _DISCARD     = {'merge_cloudnumber', 'split_cloudnumber'}
    _SUM_MERGERS = {'merge_ccs_area', 'split_ccs_area'}

    # Resolve which tracks to keep
    if keep_indices is None:
        keep_indices = np.arange(ds.sizes['tracks'])
    n_tracks = len(keep_indices)

    # track_duration for the kept tracks: read full 1D array then numpy-index
    # (1D arrays are tiny and never have the chunk-read problem)
    track_duration = ds['track_duration'].values[keep_indices].astype(int)

    # local_track_idx : 0-based index within the kept-tracks array
    # orig_track_idx  : index into the ORIGINAL ds tracks dimension
    # time_idx        : time step index within each track
    local_track_idx = np.repeat(np.arange(n_tracks), track_duration)
    orig_track_idx  = keep_indices[local_track_idx]   # map local → original
    time_idx        = np.concatenate([np.arange(d) for d in track_duration])
    n_rows          = len(orig_track_idx)

    # Maximum time step actually needed — used to clip 3D reads so we never
    # read padding time steps that belong only to excluded (e.g. tropical)
    # tracks.  time_idx.max() == track_duration.max() - 1, so this is safe.
    max_time = int(track_duration.max())
    times_total = ds.sizes.get('times', max_time)
    if max_time < times_total:
        log.info(f"  times clip: {times_total} → {max_time} "
                 f"(saves {(times_total - max_time) / times_total:.0%} on 3-D reads)")

    data = {'tracks': local_track_idx, 'relative_step': time_idx}

    # ------------------------------------------------------------------
    # Read each variable from the ORIGINAL (non-isel'd) dataset.
    #
    # ds[v].values on the original dataset triggers a single contiguous
    # netCDF4 read of the full array — no chunk-by-chunk penalty regardless
    # of whether the file is compressed/chunked.  The selection of kept
    # tracks and valid time steps is done entirely in numpy after the read.
    # ------------------------------------------------------------------
    var_times = {}
    for v in ds.data_vars:
        dims = tuple(ds[v].dims)
        t_v = time.time()

        if dims == ('tracks',):
            data[v] = ds[v].values[orig_track_idx]

        elif dims == ('tracks', 'times'):
            data[v] = ds[v].values[orig_track_idx, time_idx]

        elif dims == ('tracks', 'times', 'mergers'):
            if v in _DISCARD:
                continue                           # never read from disk
            if v in _SUM_MERGERS:
                # Clip times to max valid step — single contiguous hyperslab
                # read of (all_tracks, max_time, mergers); fill→0, sum mergers
                # Use nan_to_num before clip to also handle NaN fill values
                arr = np.nan_to_num(
                    ds[v].values[:, :max_time, :].astype(float),
                    nan=0.0, posinf=0.0, neginf=0.0,
                )
                arr = np.clip(arr, 0, None)        # numeric fill (-9999) → 0
                data[v] = arr.sum(axis=2)[orig_track_idx, time_idx]

        elif dims == ('tracks', 'times', 'nmaxpf'):
            # Clip times and take index 0 — single hyperslab read
            data[v] = ds[v].values[:, :max_time, 0][orig_track_idx, time_idx]

        # else: skip unrecognised dimensions
        var_times[v] = time.time() - t_v

    # Report slowest variables for diagnostics
    slowest = sorted(var_times.items(), key=lambda x: x[1], reverse=True)[:10]
    log.info(f"  Top-10 slowest var reads (s): "
             + ", ".join(f"{v}={t:.2f}" for v, t in slowest))

    t_df = time.time()
    df = pd.DataFrame(data)
    log.info(f"  pd.DataFrame() took {time.time() - t_df:.1f} s  "
             f"| {df.shape[0]:,} rows × {df.shape[1]} cols")
    df['relative_time_h']  = df['relative_step'] * time_res_h
    df['track_duration_h'] = df['track_duration'] * time_res_h
    df['dataset'] = name

    log.info(f"  {name:35s}: {n_tracks:>6,} tracks  →  {n_rows:>9,} valid steps "
             f"(mean {n_rows / n_tracks:.1f} h/track)")
    return df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine MCS track stats + COF 2D data for multiple sources "
                    "and save as a parquet file."
    )
    parser.add_argument('--sources', nargs='+', default=None,
                        help='Source names to process (default: all 6 sources). '
                             'Example: --sources scream IMERGv7')
    parser.add_argument('--mcs-dir', default=MCS_DIR,
                        help=f'Root directory for MCS track stats (default: {MCS_DIR})')
    parser.add_argument('--cof-dir', default=COF_DIR,
                        help=f'Directory for COF 2D netCDF files (default: {COF_DIR})')
    parser.add_argument('--output-dir', default=OUTPUT_DIR,
                        help=f'Output directory for parquet file (default: {OUTPUT_DIR})')
    parser.add_argument('--output-file', default=None,
                        help='Output parquet filename (default: '
                             'mcs_cof_trackstats_allsources.parquet)')
    parser.add_argument('--time-res-h', type=float, default=1.0,
                        help='MCS track time resolution in hours (default: 1.0)')
    parser.add_argument('--tropics-lat-threshold', type=float, default=20.0,
                        help='Absolute latitude threshold (degrees) for tropical '
                             'track removal.  Tracks with lifetime-median '
                             '|meanlat| < threshold are excluded (default: 20.0). '
                             'Set to 0 to disable.')
    return parser.parse_args()


def main():
    setup_logging()
    log = logging.getLogger(__name__)
    args = parse_args()

    sources = args.sources or DEFAULT_SOURCES
    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)

    out_file = args.output_file or "mcs_cof_trackstats_allsources.parquet"
    out_path = os.path.join(out_dir, out_file)

    log.info("=" * 70)
    log.info("MCS + COF TRACKSTATS COMBINATION  |  multi-source")
    log.info("=" * 70)
    log.info(f"Sources             : {sources}")
    log.info(f"MCS dir             : {args.mcs_dir}")
    log.info(f"COF dir             : {args.cof_dir}")
    log.info(f"Output              : {out_path}")
    log.info(f"Tropics threshold   : {args.tropics_lat_threshold}° "
             f"({'disabled' if args.tropics_lat_threshold == 0 else 'enabled'})")
    log.info("=" * 70)
    
    # ------------------------------------------------------------------
    # Step 1: Load and merge datasets
    # ------------------------------------------------------------------
    log.info("Step 1: Loading and merging track stats + COF 2D datasets...")
    ds_all = {}
    for src in sources:
        if src not in TRACKSTATS_FILES:
            log.warning(f"  {src}: no trackstats filename defined — skipping.")
            continue
        try:
            ds_all[src] = load_combined(
                src, args.mcs_dir, args.cof_dir, TRACKSTATS_FILES, logger=log
            )
        except FileNotFoundError as e:
            log.error(str(e))
            log.warning(f"  {src}: skipping due to missing file.")

    # ------------------------------------------------------------------
    # Step 1b: Compute tropical track filter indices
    # Returns integer keep_indices into the original dataset rather than
    # an isel'd dataset — this avoids lazy fancy-indexing on compressed
    # netCDF variables (which causes chunk-by-chunk reads = very slow).
    # ------------------------------------------------------------------
    keep_all = {}   # src -> np.ndarray of int (indices of kept tracks)
    if args.tropics_lat_threshold > 0:
        log.info(f"Step 1b: Filtering tropical tracks "
                 f"(|meanlat_median| < {args.tropics_lat_threshold}°)...")
        for src in list(ds_all.keys()):
            keep_all[src] = filter_extratropical_tracks(
                ds_all[src],
                tropics_lat_threshold=args.tropics_lat_threshold,
                logger=log,
            )
    else:
        log.info("Step 1b: Tropical track filtering disabled (--tropics-lat-threshold 0).")
        for src in ds_all:
            keep_all[src] = None   # None = keep all tracks

    # ------------------------------------------------------------------
    # Step 2: Convert to tidy DataFrames
    # Passes the ORIGINAL (non-isel'd) dataset to ds_to_dataframe so that
    # ds[v].values triggers a single contiguous netCDF4 read per variable.
    # The keep_indices filter is applied in numpy after the full read.
    # ------------------------------------------------------------------
    log.info(f"Step 2: Converting {len(ds_all)} sources to tidy DataFrames...")
    t0 = time.time()
    dfs = []
    for src, ds in ds_all.items():
        label = SOURCE_DISPLAY_NAMES.get(src, src)
        log.info(f"  Converting {src} ({label})...")
        t_src = time.time()
        df_src = ds_to_dataframe(ds, keep_indices=keep_all[src],
                                 time_res_h=args.time_res_h,
                                 name=label, logger=log)
        elapsed_src = time.time() - t_src
        log.info(f"  {src}: done in {elapsed_src:.1f} s "
                 f"| {df_src.shape[0]:,} rows, "
                 f"{df_src.memory_usage(deep=True).sum() / 1e6:.0f} MB")
        dfs.append(df_src)
        ds.close()

    df_all_src = pd.concat(dfs, ignore_index=True)
    elapsed = time.time() - t0

    log.info(f"Conversion time  : {elapsed:.1f} s")
    log.info(f"Combined shape   : {df_all_src.shape[0]:,} rows × {df_all_src.shape[1]} cols")
    log.info(f"Memory (approx)  : {df_all_src.memory_usage(deep=True).sum() / 1e9:.2f} GB")
    log.info(f"Sources present  : {df_all_src['dataset'].unique().tolist()}")

    # ------------------------------------------------------------------
    # Step 3: Save to parquet with snappy (fast, minimal) compression
    # ------------------------------------------------------------------
    log.info(f"Step 3: Saving combined DataFrame to parquet...")
    df_all_src.to_parquet(out_path, index=False, engine='pyarrow',
                          compression='snappy')
    size_mb = os.path.getsize(out_path) / 1e6
    log.info(f"Saved → {out_path}  ({size_mb:.1f} MB)")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    log.info("=" * 70)
    log.info("Track duration (h) summary per source:")
    df_tracks = df_all_src.drop_duplicates(subset=['dataset', 'tracks'])
    stats = (
        df_tracks.groupby('dataset')['track_duration_h']
        .describe(percentiles=[0.5, 0.90, 0.95])
        [['count', '50%', '90%', '95%', 'max']]
    )
    for line in stats.to_string().split('\n'):
        log.info(f"  {line}")
    log.info("=" * 70)
    log.info("DONE")


if __name__ == "__main__":
    main()
