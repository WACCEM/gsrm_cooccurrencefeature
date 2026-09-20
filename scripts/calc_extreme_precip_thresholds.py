"""
Calculate extreme precipitation percentiles from HEALPix grid data.

This script computes precipitation percentiles (e.g., 90th, 95th) at different
time scales (6-hourly, daily, etc.) from Zarr datasets on HEALPix grid.

By default each source reads its own 6-hourly precipitation product. With --input_zarr the input is any Zarr store on the
HEALPix grid, for example Step 1's tot_pr (the 6-hourly window-mean precipitation the attribution uses, so thresholds and attributed
precipitation come from one field):

  --input_zarr <data root>/mcs_masks/<source>_mcs_masks_hp8.zarr --input_var tot_pr

Author: Zhe Feng, zhe.feng@pnnl.gov
"""
import numpy as np
import sys
import os
from pathlib import Path
import yaml
import xarray as xr
import pandas as pd
import time
import argparse
import intake
import logging
import math
import easygems.healpix as egh
sys.path.append(str(Path(__file__).parent.parent))
from src.cof_paths import data_root


def parse_cmd_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description='Calculate extreme precipitation percentiles')
    parser.add_argument('--catalog_source', type=str, required=True,
                        help='Catalog source name (e.g., scream_ne120, IR_IMERG)')
    parser.add_argument('--config_file', type=str, default='/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources.yaml',
                        help='Path to configuration YAML file (default: /global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources.yaml)')
    parser.add_argument('--zoom', type=int, default=8,
                        help='HEALPix zoom level (default: 8)')
    parser.add_argument('--percentiles', type=float, nargs='+', default=[90, 95],
                        help='Percentiles to compute (default: 90 95)')
    parser.add_argument('--time_durations', type=str, nargs='+', default=['6h'],
                        help='Time durations for resampling (default: 6h)')
    parser.add_argument('--method', type=str, default='linear',
                        choices=['linear', 'lower', 'higher', 'midpoint', 'nearest'],
                        help='Quantile interpolation method (default: linear)')
    parser.add_argument('--min_precip_threshold', type=float, default=0.01,
                        help='Minimum precipitation threshold (mm/h) to exclude from percentile calculation. '
                             'Values below this threshold are treated as missing. (default: 0.01 mm/h)')
    parser.add_argument('--input_zarr', type=str, default=None,
                        help='Explicit path to a Zarr store on the HEALPix grid that replaces the source\'s normal input (catalog '
                             'or 6-hourly product), for ANY source. The store\'s zoom must match --zoom (checked). All of its '
                             'time steps are used unless --start_time/--end_time are given (the config dates apply only to the '
                             'normal input). Example: Step 1\'s <data root>/mcs_masks/<source>_mcs_masks_hp8.zarr with '
                             '--input_var tot_pr, or the non-IR IMERG 6-hourly store, which has data at all latitudes.')
    parser.add_argument('--input_var', type=str, default=None,
                        help='Variable to read from --input_zarr (default: the config\'s varname_precip_liq). When given, the '
                             'field is used as it is (no frozen precipitation is added) and --input_factor defaults to 1.')
    parser.add_argument('--input_factor', type=float, default=None,
                        help='Factor that converts the input variable to mm/h (default: the config\'s pr_convert_factor, or 1 '
                             'when --input_var is given, e.g. tot_pr is already in mm/h)')
    parser.add_argument('--start_time', type=str, default=None,
                        help='Inclusive start of the period to use, any label xarray accepts for .sel(time=slice(...)), e.g. '
                             '"2019-08-01T00" or "2019-08". Default: the config start date for the normal input, the whole '
                             'store for --input_zarr.')
    parser.add_argument('--end_time', type=str, default=None,
                        help='Inclusive end of the period to use; a partial string such as "2020-08" extends to the end of that '
                             'month. Default: the config end date for the normal input, the whole store for --input_zarr.')
    parser.add_argument('--cell_chunk_size', type=int, default=None,
                        help='Number of HEALPix cells processed per block when computing the quantiles over time (all time steps '
                             'x this many cells are held in memory at once). Default: sized from the available memory and the '
                             'store\'s native chunk width. Pass a value to cap the memory when other jobs share the node.')
    parser.add_argument('--available_memory_gb', type=float, default=None,
                        help='Memory (GB) to budget for the block size when --cell_chunk_size is not given (default: the '
                             'available memory of the node, detected with psutil)')
    parser.add_argument('--output_dir', type=str, default=None,
                        help='Output directory for NetCDF files (default: extreme_precip/ under the pipeline data root, see '
                             'src/cof_paths.py; production unless COF_DATA_ROOT is set)')
    parser.add_argument('--version', type=str, default='v1',
                        help='Version string for output files (default: v1)')
    
    args = parser.parse_args()
    args_dict = vars(args)
    return args_dict


def setup_logging():
    """
    Set up logging configuration.
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)


def load_config(config_file, catalog_source):
    """
    Load configuration from YAML file for a specific catalog source.
    
    Args:
        config_file: str
            Path to the YAML configuration file
        catalog_source: str
            The catalog source key to load configuration for
            
    Returns:
        dict: Configuration dictionary for the specified source
    """
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    if catalog_source not in config:
        raise ValueError(f"Catalog source '{catalog_source}' not found in config file. "
                        f"Available sources: {list(config.keys())}")
    
    return config[catalog_source]


#--------------------------------------------------------------------------------------------------
# The helpers below are the ones of calc_extreme_precip_thresholds_1h.py (copied unchanged): the time subset, the zoom check of an
# explicit --input_zarr, and the memory-bounded quantile over time that processes the cells in blocks.
#--------------------------------------------------------------------------------------------------

def _time_to_isoformat(t):
    """
    Convert a single time coordinate value to an ISO-8601 string, regardless of
    whether it uses a standard (numpy.datetime64) or a non-standard-calendar
    (cftime, e.g. SCREAM's 'noleap'/365_day) representation.

    pd.Timestamp() cannot represent non-standard calendars at all and raises
    TypeError on cftime objects; those objects already provide their own working
    isoformat(), so this only falls back to pd.Timestamp() for values that lack
    one (e.g. numpy.datetime64, which has no isoformat() method of its own).
    """
    if hasattr(t, 'isoformat'):
        return t.isoformat()
    return pd.Timestamp(t).isoformat()


def subset_time(ds_p, start_time=None, end_time=None, logger=None):
    """
    Restrict a dataset to an inclusive time period, if requested.

    Args:
        ds_p: xr.Dataset
            Source dataset, as returned by the open step in load_precipitation_data()
        start_time: str, optional
            Inclusive start bound (any label xarray accepts for time slicing, e.g.
            '2018-01-01T00', or a partial string like '2018'). None: no lower bound.
        end_time: str, optional
            Inclusive end bound. A partial string extends to the end of the implied
            interval (e.g. '2022' includes all of 2022), matching xarray's own
            partial-string slicing semantics. None: no upper bound.
        logger: logging.Logger
            Logger instance

    Returns:
        xr.Dataset: ds_p unchanged if both bounds are None or 'time' is not a
            dimension; otherwise ds_p.sel(time=slice(start_time, end_time))
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    if start_time is None and end_time is None:
        return ds_p

    if 'time' not in ds_p.dims:
        logger.warning("--start_time/--end_time given but dataset has no 'time' "
                        "dimension; ignoring.")
        return ds_p

    n_before = ds_p.sizes['time']
    logger.info(f"Subsetting time period: start_time={start_time}, end_time={end_time}")
    ds_sub = ds_p.sel(time=slice(start_time, end_time))
    n_after = ds_sub.sizes['time']

    if n_after == 0:
        avail_start = _time_to_isoformat(ds_p['time'].values[0])
        avail_end = _time_to_isoformat(ds_p['time'].values[-1])
        raise ValueError(f"Requested time period [{start_time}, {end_time}] does not "
                         f"overlap the dataset's available range "
                         f"[{avail_start}, {avail_end}].")

    logger.info(f"Time steps: {n_before} -> {n_after} "
                f"({_time_to_isoformat(ds_sub['time'].values[0])} to "
                f"{_time_to_isoformat(ds_sub['time'].values[-1])})")

    return ds_sub


def _validate_zoom(ds_p, zoom, source_desc, logger=None):
    """
    Confirm an explicitly-loaded Dataset's actual HEALPix zoom matches --zoom.

    Only relevant for --input_zarr, where data no longer comes from --zoom (unlike
    the catalog/IMERG/GSMAP paths, which select data BY zoom and so are consistent
    by construction). --zoom still names the output file
    (f"..._hp{zoom}_{version}.nc") and its zoom_level global attribute, so a
    mismatch would silently mislabel the output -- or overwrite a good file at
    another zoom -- rather than failing. HEALPix zoom and nside are related by
    nside = 2**zoom (confirmed for this repo's stores: e.g. casesm2_10km_nocumulus's
    healpix_nside=512 at zoom=9 -> log2(512)=9).

    Args:
        ds_p: xr.Dataset
            Dataset just opened from the explicit store (before or after
            egh.attach_coords() -- both leave the 'crs' variable's
            healpix_nside attribute, which this reads, untouched)
        zoom: int
            The --zoom the caller requested
        source_desc: str
            Path or other description of the store, for the error message
        logger: logging.Logger
            Logger instance

    Raises:
        ValueError: if the store's zoom cannot be determined, or does not match
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    try:
        nside = egh.get_nside(ds_p)
        actual_zoom = int(round(math.log2(nside)))
    except Exception as exc:
        raise ValueError(f"Could not determine the HEALPix zoom of --input_zarr "
                        f"'{source_desc}' to validate against --zoom {zoom}: {exc}")

    if actual_zoom != zoom:
        raise ValueError(f"--input_zarr '{source_desc}' is at HEALPix zoom "
                        f"{actual_zoom} (nside={nside}), but --zoom {zoom} was "
                        f"requested. --zoom names the output file and its "
                        f"zoom_level attribute, so these must match -- pass "
                        f"--zoom {actual_zoom} instead.")
    logger.info(f"Confirmed --input_zarr matches requested zoom: {actual_zoom}")


def determine_cell_chunk_size(pr_resampled, cell_chunk_size=None, available_memory_gb=None,
                              mem_safety_fraction=0.6, overhead_factor=2.5, min_chunk_size=1024,
                              logger=None):
    """
    Determine how many HEALPix cells to process per block in _quantile_over_time_blocked(),
    sized to (a) fit safely in available memory and (b) align to the dataset's native
    Zarr chunk boundaries to minimize redundant chunk decompression.

    Why alignment matters: a compressed Zarr chunk must be decompressed whole (blosc
    has no partial/random-access decompression within a chunk). If a block only needs
    a fraction of a native chunk's cells, every *other* block sharing that same native
    chunk re-decompresses it independently from scratch -- e.g. a block 1/16th the
    width of the native chunk causes that native chunk to be decompressed 16 times
    total (once per sibling block) to serve the full cell range, even though the
    underlying bytes only need decompressing once. This is the dominant cost, not raw
    I/O bandwidth: a too-small cell_chunk_size (e.g. a small fixed default) makes each
    block correct and memory-safe, but wastefully slow. Sizing a block to a whole
    native chunk (or a whole multiple of it) eliminates the redundancy entirely;
    sizing it to an exact divisor of the native chunk minimizes -- but can't
    eliminate -- the redundancy when memory doesn't allow a full native chunk per
    block.

    Args:
        pr_resampled: xr.DataArray
            The (already prepared) precipitation array this sizing is for. Its
            actual 'time' length and 'cell' dask chunking are read directly here
            (not hardcoded), so the same call automatically adapts to HEALPix zoom
            level (more/fewer total cells, and possibly a different native chunk
            width) and to the time period being processed -- e.g. a single year in
            the annual loop naturally gets a far larger safe chunk size than the
            multi-decade all-period call, because n_time is much smaller.
        cell_chunk_size: int, optional
            If given, used as-is (manual override, no auto-tuning). If None
            (default), computed dynamically as described above.
        available_memory_gb: float, optional
            System memory to budget against, in GB. If None, auto-detected via
            psutil.virtual_memory().available. Pass explicitly when running under a
            job scheduler (e.g. SLURM) if the auto-detected node-level value doesn't
            reliably reflect the job's actual memory allocation/cgroup limit.
        mem_safety_fraction: float
            Fraction of available_memory_gb to budget for one block. Blocks are
            processed strictly one at a time (no concurrency), so this only needs to
            leave headroom for the Python/dask/OS baseline and other steps in the
            script, not for multiple blocks in flight. Default 0.6.
        overhead_factor: float
            Multiplier from a block's raw byte size (n_time * cell_chunk_size *
            dtype size) to its actual peak memory while being processed (the
            threshold-masking copy plus quantile's internal sort scratch). Default
            2.5, based on empirical observation.
        min_chunk_size: int
            Floor on the returned chunk size regardless of how tight the memory
            budget is. Default 1024.
        logger: logging.Logger
            Logger instance

    Returns:
        int: number of cells to process per block
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    ncells = pr_resampled.sizes['cell']

    if cell_chunk_size is not None:
        logger.info(f"Using explicit cell_chunk_size={cell_chunk_size} (auto-sizing skipped)")
        return min(cell_chunk_size, ncells)

    n_time = pr_resampled.sizes['time']
    dtype_bytes = pr_resampled.dtype.itemsize
    cell_axis = pr_resampled.dims.index('cell')
    # Native Zarr chunk width along 'cell', read directly from the actual dask
    # chunking rather than assumed -- correct regardless of HEALPix zoom level or
    # future changes to how the source Zarr stores are chunked.
    native_cell_chunk = pr_resampled.chunks[cell_axis][0] if pr_resampled.chunks is not None else ncells

    if available_memory_gb is not None:
        available_bytes = available_memory_gb * 1024**3
    else:
        try:
            import psutil
            available_bytes = psutil.virtual_memory().available
            logger.info(f"Auto-detected available memory: {available_bytes / 1e9:.1f} GB "
                        f"(pass --available_memory_gb to override, e.g. if running under a "
                        f"SLURM allocation where node-level auto-detection may be unreliable)")
        except ImportError:
            available_bytes = 64 * 1024**3
            logger.warning("psutil not installed; assuming 64GB available memory for automatic "
                            "chunk sizing. Install psutil, or pass --available_memory_gb / "
                            "--cell_chunk_size explicitly, for accurate sizing.")

    raw_budget_bytes = (available_bytes * mem_safety_fraction) / overhead_factor
    max_cells_by_memory = max(int(raw_budget_bytes / (n_time * dtype_bytes)), min_chunk_size)

    if max_cells_by_memory >= native_cell_chunk:
        # Enough memory for one or more FULL native chunks per block -- round down
        # to the largest whole multiple of the native chunk width (zero redundant reads)
        n_native_chunks_per_block = max(1, max_cells_by_memory // native_cell_chunk)
        chunk_size = n_native_chunks_per_block * native_cell_chunk
        redundancy = 1
    else:
        # Not enough memory for one full native chunk -- pick the largest divisor of
        # the native chunk width that still fits, keeping each block's reads
        # confined to a single native chunk (minimizes, though doesn't eliminate,
        # cross-block redundant decompression)
        divisor = native_cell_chunk
        while divisor > max_cells_by_memory and divisor > min_chunk_size:
            divisor //= 2
        chunk_size = max(divisor, min_chunk_size)
        redundancy = max(1, native_cell_chunk // chunk_size)

    chunk_size = min(chunk_size, ncells)
    n_blocks = -(-ncells // chunk_size)  # ceil division
    est_raw_gb = n_time * chunk_size * dtype_bytes / 1e9
    logger.info(f"Auto-sized cell_chunk_size={chunk_size} ({n_blocks} block(s)) for "
                f"n_time={n_time}, native_cell_chunk={native_cell_chunk}: "
                f"~{est_raw_gb:.1f}GB raw (~{est_raw_gb * overhead_factor:.1f}GB peak) per "
                f"block, {redundancy}x redundant native-chunk re-decompression")

    return chunk_size


def _quantile_over_time_blocked(pr_resampled, percentiles, method='linear',
                                cell_chunk_size=None, available_memory_gb=None, logger=None):
    """
    Compute quantile(s) over 'time' by looping over 'cell' in fixed-size blocks and
    materializing one block at a time. Numerically identical to
    pr_resampled.quantile([p/100 for p in percentiles], dim='time', method=method,
    skipna=True) (differs only by float32 round-off), but exists purely to bound
    peak memory.

    Why this is necessary: xr.DataArray.quantile() delegates to
    dask.array.nanquantile, which -- whenever 'time' spans more than one chunk
    (always true for these multi-year Zarr stores, natively chunked ~24 steps at a
    time) -- rechunks 'time' into a single chunk and picks an 'auto' size for the
    other dimension targeting dask's default 128MiB/chunk. For a multi-decade hourly
    record that produces thousands of tiny output blocks, each requiring the full
    time axis to be re-read from every native chunk (e.g. 24 years of IMERG at zoom
    8 explodes to ~4,900 output blocks and a ~274,000-task graph), which is what
    exhausts memory even on a ~500GB node. Processing 'cell' in controlled,
    sequential blocks instead keeps peak memory bounded to roughly one block's own
    footprint and sidesteps the pathological fan-out entirely. See
    determine_cell_chunk_size() for how the block size itself is chosen -- it is
    sized dynamically by default to also minimize the redundant chunk
    re-decompression that a too-small fixed block size would otherwise cause.

    Args:
        pr_resampled: xr.DataArray
            Precipitation data array with dimensions (time, cell), already prepared
            (threshold applied, resampled) via prepare_precip()
        percentiles: list of float
            Percentiles to compute (e.g., [90, 95] for 90th and 95th percentiles)
        method: str
            Interpolation method for quantile calculation. Options include:
            - 'linear': linear interpolation (default)
            - 'lower': lower value
            - 'higher': higher value
            - 'midpoint': midpoint of two nearest values
            - 'nearest': nearest value
        cell_chunk_size: int, optional
            Number of cells to process per block. If None (default), determined
            dynamically per-call via determine_cell_chunk_size() -- see there for
            details. Pass an explicit value to override.
        available_memory_gb: float, optional
            Forwarded to determine_cell_chunk_size() when cell_chunk_size is None.
        logger: logging.Logger
            Logger instance

    Returns:
        dict: Dictionary with percentile values as keys and already-computed
              (eager) xr.DataArrays (dim: cell) as values, e.g., {90: pr_p90, 95: pr_p95}
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    cell_chunk_size = determine_cell_chunk_size(pr_resampled, cell_chunk_size=cell_chunk_size,
                                                available_memory_gb=available_memory_gb, logger=logger)

    ncells = pr_resampled.sizes['cell']
    cell_coord = pr_resampled['cell']
    q_values = [p / 100.0 for p in percentiles]

    out_arrays = {p: np.full(ncells, np.nan, dtype='float32') for p in percentiles}

    n_blocks = -(-ncells // cell_chunk_size)  # ceil division
    for i, start in enumerate(range(0, ncells, cell_chunk_size)):
        end = min(start + cell_chunk_size, ncells)
        logger.info(f"  Block {i + 1}/{n_blocks}: cells {start}-{end - 1} of {ncells}...")
        # A single, bounded .load() -- materializes just this cell block, all time steps
        block = pr_resampled.isel(cell=slice(start, end)).load()
        # np.nanquantile skips NaN values, matching the skipna=True behavior used
        # elsewhere in this script. Computing all requested percentiles together
        # avoids reloading the block once per percentile.
        block_result = np.nanquantile(block.values, q_values, axis=0, method=method)
        for qi, p in enumerate(percentiles):
            out_arrays[p][start:end] = block_result[qi]
        del block, block_result

    results = {}
    for p in percentiles:
        results[p] = xr.DataArray(out_arrays[p], dims=['cell'], coords={'cell': cell_coord})

    return results


def calc_precip_percentiles(pr, percentiles=[90, 95], time_duration='6h', method='linear', 
                           min_precip_threshold=None, logger=None, blocked=True,
                           cell_chunk_size=None, available_memory_gb=None):
    """
    Calculate precipitation percentiles at specified time scales.
    
    Args:
        pr: xr.DataArray
            Precipitation data array with dimensions (time, cell) in mm/h
        percentiles: list of float
            Percentiles to compute (e.g., [90, 95] for 90th and 95th percentiles)
        time_duration: str
            Time duration for resampling. Options include:
            - '6h' or '6H': 6-hourly (no resampling if already 6-hourly)
            - '1D' or '24h': daily
            - '12h' or '12H': 12-hourly
            Default is '6h' (keeps original 6-hourly resolution)
        method: str
            Interpolation method for quantile calculation. Options include:
            - 'linear': linear interpolation (default)
            - 'lower': lower value
            - 'higher': higher value
            - 'midpoint': midpoint of two nearest values
            - 'nearest': nearest value
        min_precip_threshold: float, optional
            Minimum precipitation threshold (mm/h). Values below this threshold
            are excluded from percentile calculations (treated as NaN).
            If None, all precipitation values are included. Default is None.
        logger: logging.Logger
            Logger instance
        blocked: bool
            True (default): compute the quantiles over time in blocks of cells, which bounds the memory (see
            _quantile_over_time_blocked); the values equal those of xarray's quantile up to float32 round-off.
            False: xarray's quantile on the whole array (the previous implementation, kept for comparison).
        cell_chunk_size: int, optional
            Cells per block when blocked; default: sized from the available memory and the native chunk width.
        available_memory_gb: float, optional
            Memory to budget for the block size when cell_chunk_size is not given.
            
    Returns:
        dict: Dictionary with percentile values as keys and DataArrays as values
              e.g., {90: pr_p90, 95: pr_p95}
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Apply minimum precipitation threshold if specified
    if min_precip_threshold is not None:
        logger.info(f"Applying minimum precipitation threshold: {min_precip_threshold} mm/h")
        logger.info(f"Values below threshold will be excluded from percentile calculation")
        # Replace values below threshold with NaN
        pr_filtered = pr.where(pr >= min_precip_threshold)
        # # Count excluded values
        # n_below = int((pr < min_precip_threshold).sum().values)
        # n_total = int(pr.size)
        # pct_excluded = 100.0 * n_below / n_total if n_total > 0 else 0
        # logger.info(f"Excluded {n_below:,} values ({pct_excluded:.2f}%) below threshold")
    else:
        logger.info("No minimum precipitation threshold applied")
        pr_filtered = pr
    
    # Resample to specified time duration if needed
    if time_duration.lower() not in ['6h', '6hr']:
        # Convert time duration to pandas-compatible frequency string
        freq = time_duration.upper().replace('HR', 'H')
        logger.info(f"Resampling precipitation from 6-hourly to {freq}...")
        # Resample by taking mean precipitation rate (ignoring NaN)
        pr_resampled = pr_filtered.resample(time=freq).mean(keep_attrs=True)
    else:
        logger.info(f"Using original 6-hourly precipitation data...")
        pr_resampled = pr_filtered
    
    # Calculate percentiles. Blocked: the cells are processed in blocks, all percentiles together, NaN values ignored.
    results = {}
    blocked_results = None
    if blocked:
        logger.info(f"Computing the {percentiles} percentiles over time in blocks of cells...")
        blocked_results = _quantile_over_time_blocked(pr_resampled, percentiles, method=method,
                                                      cell_chunk_size=cell_chunk_size,
                                                      available_memory_gb=available_memory_gb, logger=logger)
    for p in percentiles:
        if blocked:
            pr_percentile = blocked_results[p]
        else:
            logger.info(f"Computing {p}th percentile...")
            quantile_value = p / 100.0
            # skipna=True by default in quantile(), so NaN values are ignored
            pr_percentile = pr_resampled.quantile(quantile_value, dim='time', method=method, skipna=True)
        pr_percentile.name = f'pr_p{p}'
        pr_percentile.attrs['long_name'] = f'{p}th percentile precipitation'
        pr_percentile.attrs['units'] = 'mm/h'
        pr_percentile.attrs['time_duration'] = time_duration
        pr_percentile.attrs['method'] = method
        if min_precip_threshold is not None:
            pr_percentile.attrs['min_precip_threshold'] = f'{min_precip_threshold} mm/h'
            pr_percentile.attrs['note'] = f'Precipitation below {min_precip_threshold} mm/h excluded from calculation'
        results[p] = pr_percentile
    
    return results


def write_netcdf(results_dict, ds_p, output_filename, zoom, source_name, 
                 start_datetime, end_datetime, time_duration, method, 
                 min_precip_threshold=None, logger=None, extra_attrs=None):
    """
    Write precipitation percentiles to a NetCDF file.
    
    Args:
        results_dict: dict
            Dictionary with percentile values as keys and DataArrays as values
        ds_p: xr.Dataset
            Original precipitation dataset (for coordinates)
        output_filename: str
            Output NetCDF filename
        zoom: int
            HEALPix zoom level
        source_name: str
            Data source name
        start_datetime: str
            Start date/time of data
        end_datetime: str
            End date/time of data
        time_duration: str
            Time duration used for resampling
        method: str
            Quantile interpolation method used
        min_precip_threshold: float, optional
            Minimum precipitation threshold used
        logger: logging.Logger
            Logger instance
        extra_attrs: dict, optional
            Further global attributes to record, e.g. where the precipitation came from
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    logger.info(f'Preparing data for output file: {output_filename}')
    
    # Create variables dictionary
    var_dict = {}
    for p, data in results_dict.items():
        var_name = f'pr_p{int(p)}'
        var_dict[var_name] = (['cell'], data.values)
    
    # Create coordinates
    coord_dict = {
        'cell': (['cell'], ds_p['cell'].values),
        'lat': (['cell'], ds_p['lat'].values),
        'lon': (['cell'], ds_p['lon'].values),
    }
    
    # Add crs if available
    if 'crs' in ds_p:
        coord_dict['crs'] = ds_p['crs'].values
    
    # Create global attributes
    gattr_dict = {
        'Title': 'Precipitation percentiles at different time scales',
        'contact': 'Zhe Feng, zhe.feng@pnnl.gov',
        'source_name': source_name,
        'start_date': start_datetime,
        'end_date': end_datetime,
        'created_on': time.ctime(time.time()),
        'grid_type': 'HEALPix',
        'zoom_level': zoom,
        'time_duration': time_duration,
        'quantile_method': method,
    }
    
    if min_precip_threshold is not None:
        gattr_dict['min_precip_threshold'] = f'{min_precip_threshold} mm/h'
        gattr_dict['threshold_note'] = f'Precipitation below {min_precip_threshold} mm/h excluded from percentile calculation'
    if extra_attrs:
        gattr_dict.update(extra_attrs)
    
    # Create output dataset
    dsout = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)
    
    # Add coordinate attributes
    dsout['cell'].attrs['long_name'] = 'HEALPix cell index'
    dsout['lon'].attrs['long_name'] = 'Longitude'
    dsout['lon'].attrs['units'] = 'degree'
    dsout['lat'].attrs['long_name'] = 'Latitude'
    dsout['lat'].attrs['units'] = 'degree'
    
    # Add variable attributes
    for p in results_dict.keys():
        var_name = f'pr_p{int(p)}'
        dsout[var_name].attrs['long_name'] = f'{int(p)}th percentile precipitation'
        dsout[var_name].attrs['units'] = 'mm/h'
        dsout[var_name].attrs['time_duration'] = time_duration
        dsout[var_name].attrs['percentile'] = p
        dsout[var_name].attrs['method'] = method
        if min_precip_threshold is not None:
            dsout[var_name].attrs['min_precip_threshold'] = f'{min_precip_threshold} mm/h'
    
    # Save the output file
    fillvalue = np.nan
    comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}
    
    logger.info(f'Writing output file: {output_filename}')
    dsout.to_netcdf(path=output_filename, mode='w', format='NETCDF4', encoding=encoding)
    logger.info(f'Successfully wrote: {output_filename}')
    
    return dsout


def _check_precip_scale(pr, logger, low=1.0e-3, high=20.0):
    """
    Guard against a wrong --input_var or --input_factor: the domain mean of the precipitation, over a few time steps spread over the
    record, must be a plausible rate in mm/h (typically about 0.1). Raises ValueError otherwise, e.g. for a field still in m/s or
    mm/day, before a wrong threshold file is written.
    """
    n_time = pr.sizes['time']
    idx = np.unique(np.linspace(0, n_time - 1, min(n_time, 8)).astype(int))
    means = [float(pr.isel(time=int(i)).mean(skipna=True)) for i in idx]
    means = [m for m in means if np.isfinite(m)]
    if not means:
        logger.warning("Scale check skipped: the sampled time steps are all NaN")
        return
    mean = float(np.mean(means))
    if not (low <= mean <= high):
        raise ValueError(f"The domain-mean precipitation of the input is {mean:.3g}, outside the plausible {low}-{high} mm/h. "
                         f"Check --input_var and --input_factor (the field must be in mm/h after the factor is applied).")
    logger.info(f"Input scale check passed: domain mean {mean:.4f} mm/h over {len(means)} sampled time steps")


def load_precipitation_data(config_file, catalog_source, zoom, logger=None, input_zarr=None, input_var=None, input_factor=None):
    """
    Load precipitation data from catalog or direct Zarr file.
    
    Args:
        config_file: str
            Path to configuration YAML file
        catalog_source: str
            Catalog source name
        zoom: int
            HEALPix zoom level
        logger: logging.Logger
            Logger instance
        input_zarr: str, optional
            Explicit Zarr store that replaces the source's normal input (see the module docstring). Its HEALPix zoom must
            equal zoom.
        input_var: str, optional
            Variable to read from input_zarr (default: the config's varname_precip_liq). When given, no frozen precipitation
            is added and input_factor defaults to 1.
        input_factor: float, optional
            Conversion factor to mm/h (default: the config's pr_convert_factor, or 1 when input_var is given)
            
    Returns:
        tuple: (pr DataArray, ds_p Dataset, config dict). pr is (time, cell) in mm/h; for input_zarr its attrs record
        'input_var' and 'input_factor'.
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Load configuration for the specified source
    config = load_config(config_file, catalog_source)
    
    # Extract configuration variables
    source_name = config.get('source_name')
    catalog_location = config.get('catalog_location', 'NERSC')
    catalog_params = config.get('catalog_params', {}).copy()
    varname_precip_liq = config.get('varname_precip_liq')
    varname_precip_ice = config.get('varname_precip_ice')
    pr_convert_factor = config.get('pr_convert_factor')
    start_datetime = config.get('start_datetime')
    end_datetime = config.get('end_datetime')
    
    # Catalog parameters
    catalog_file = "https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml"
    # catalog_file = "/global/homes/f/feng045/program/hackathon/catalog/NERSC/main.yaml"
    
    if input_zarr is not None:
        # An explicit store replaces the source's normal input, whatever the source is
        if not os.path.exists(input_zarr):
            raise FileNotFoundError(f"--input_zarr does not exist: {input_zarr}")
        logger.info(f"Loading {catalog_source} from explicit --input_zarr (NOT the source's normal input): {input_zarr}")
        # consolidated=None: use the consolidated metadata when present, scan the store otherwise
        ds_p = xr.open_zarr(input_zarr, consolidated=None)
        ds_p = ds_p.pipe(egh.attach_coords)
        # --zoom names the output file and its zoom_level attribute, so a store at another zoom must not pass silently
        _validate_zoom(ds_p, zoom, input_zarr, logger=logger)

    # Special treatment for certain datasets not in the catalog
    elif catalog_source == "IR_IMERG":
        # Special case for IMERG data (not in catalog yet)
        dir_healpix = "/pscratch/sd/w/wcmca1/GPM/healpix/"
        in_basename = f"IMERG_V7_"
        time_res = "6H"
        in_zarr = f"{dir_healpix}{in_basename}{time_res}_zoom{zoom}_20190101_20211231.zarr"
        # Read IMERG dataset
        logger.info(f"Loading IMERG dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)
    
    elif catalog_source == "scream_ne120":
        dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/scream/"
        in_basename = f"scream_pr"
        time_res = "6h"
        in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
        # Read SCREAM dataset
        logger.info(f"Loading SCREAM dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)
    
    elif catalog_source == "nicam_gl11":
        dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/nicam_gl11/shifted/"
        in_basename = f"NICAM_pr"
        time_res = "6h"
        in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
        # Read NICAM dataset
        logger.info(f"Loading NICAM dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)
    
    elif catalog_source == "um_glm_n2560_RAL3p3":
        dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/um_glm_n2560_RAL3p3/"
        in_basename = f"um_glm_n2560_RAL3p3_pr"
        time_res = "6h"
        in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
        # Read UM dataset
        logger.info(f"Loading UM dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)
    
    elif catalog_source == "casesm2_10km_nocumulus":
        dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/casesm2_10km_nocumulus/"
        in_basename = f"casesm2_10km_nocumulus_pr"
        time_res = "6h"
        in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
        # Read CASESM2 dataset
        logger.info(f"Loading CASESM2 dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)
    
    else:
        # Load the HEALPix catalog
        logger.info(f"Loading HEALPix catalog: {catalog_file}")
        in_catalog = intake.open_catalog(catalog_file)
        if catalog_location:
            in_catalog = in_catalog[catalog_location]
        
        # Get the DataSet from the catalog
        ds_p = in_catalog[catalog_source](**catalog_params).to_dask()
        # Add lat/lon coordinates to the HEALPix DataSet
        ds_p = ds_p.pipe(egh.attach_coords)
    
    # Which variable, which conversion factor, and whether frozen precipitation is added. Naming --input_var means the field is
    # used as it is (tot_pr is already in mm/h and already the total), so factor 1 and no frozen precipitation.
    if input_var is not None:
        varname_precip_liq = input_var
        varname_precip_ice = None
        pr_convert_factor = 1.0 if input_factor is None else input_factor
    elif input_factor is not None:
        pr_convert_factor = input_factor
    
    # Check liquid precipitation variable
    if varname_precip_liq not in list(ds_p.keys()):
        raise KeyError(f"Variable '{varname_precip_liq}' not found in the input for '{catalog_source}'. "
                       f"Available variables: {sorted(ds_p.keys())}")
    # Convert liquid precipitation to mm/h
    pr = ds_p[varname_precip_liq] * pr_convert_factor
    
    # Check if the ice precipitation variable exist in the dataset
    if varname_precip_ice in list(ds_p.keys()):
        # Convert ice precipitation to liquid equivalent
        prs = ds_p[varname_precip_ice] * pr_convert_factor
        # Add ice precipitation to get total precipitation
        pr = pr + prs
    
    # The quantiles are taken over axis 0, so time must come first
    if set(pr.dims) == {'time', 'cell'}:
        pr = pr.transpose('time', 'cell')
    
    if input_zarr is not None:
        _check_precip_scale(pr, logger)
        pr.attrs['input_var'] = varname_precip_liq
        pr.attrs['input_factor'] = float(pr_convert_factor)
    
    logger.info(f"Precipitation data loaded: {pr.shape}")
    
    return pr, ds_p, config


def main():
    """
    Main function to calculate extreme precipitation percentiles.
    """
    # Parse command-line arguments
    args_dict = parse_cmd_args()
    
    # Setup logging
    logger = setup_logging()
    logger.info("=" * 60)
    logger.info("Starting extreme precipitation percentile calculation")
    logger.info("=" * 60)
    
    # Extract arguments
    config_file = args_dict['config_file']
    catalog_source = args_dict['catalog_source']
    zoom = args_dict['zoom']
    percentiles = args_dict['percentiles']
    time_durations = args_dict['time_durations']
    method = args_dict['method']
    min_precip_threshold = args_dict['min_precip_threshold']
    output_dir = args_dict['output_dir']
    version = args_dict['version']
    input_zarr = args_dict['input_zarr']
    input_var = args_dict['input_var']
    input_factor = args_dict['input_factor']
    start_time = args_dict['start_time']
    end_time = args_dict['end_time']
    cell_chunk_size = args_dict['cell_chunk_size']
    available_memory_gb = args_dict['available_memory_gb']

    if (input_var is not None or input_factor is not None) and input_zarr is None:
        sys.exit("--input_var and --input_factor only apply together with --input_zarr")
    
    # Default output directory: extreme_precip/ under the pipeline data root (src/cof_paths.py)
    if output_dir is None:
        output_dir = f"{data_root(logger)}extreme_precip/"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load precipitation data
    logger.info(f"Loading precipitation data for {catalog_source}...")
    pr, ds_p, config = load_precipitation_data(config_file, catalog_source, zoom, logger,
                                               input_zarr=input_zarr, input_var=input_var, input_factor=input_factor)

    source_name = config.get('source_name')
    start_datetime = config.get('start_datetime')
    end_datetime = config.get('end_datetime')

    # Subset precipitation to the period asked for: --start_time/--end_time, else the configured time range of the normal
    # input; an explicit --input_zarr without times uses every time step it has
    if start_time is not None or end_time is not None:
        pr = subset_time(pr, start_time, end_time, logger)
    elif input_zarr is None and start_datetime and end_datetime:
        logger.info(f"Subsetting time range: {start_datetime} to {end_datetime}")
        pr = pr.sel(time=slice(str(start_datetime), str(end_datetime)))
        logger.info(f"Time steps after subsetting: {len(pr.time)} "
                    f"({pr.time.values[0]} to {pr.time.values[-1]})")

    extra_attrs = None
    if input_zarr is not None:
        # The period actually used goes into the file, and where the precipitation came from
        start_datetime = _time_to_isoformat(pr.time.values[0])
        end_datetime = _time_to_isoformat(pr.time.values[-1])
        logger.info(f"Input time steps: {pr.sizes['time']} ({start_datetime} to {end_datetime})")
        extra_attrs = {
            'input_zarr': input_zarr,
            'input_var': pr.attrs.get('input_var', ''),
            'input_factor': pr.attrs.get('input_factor', 1.0),
            'frames_in_input': int(pr.sizes['time']),
            'precipitation_source': f"variable {pr.attrs.get('input_var', '')} of {input_zarr} "
                                    f"({pr.sizes['time']} time steps, mm/h after the factor)",
        }

    # Process each time duration
    for time_duration in time_durations:
        logger.info("=" * 60)
        logger.info(f"Processing time duration: {time_duration}")
        logger.info("=" * 60)
        
        # Calculate percentiles
        results = calc_precip_percentiles(pr, percentiles=percentiles, 
                                         time_duration=time_duration,
                                         method=method,
                                         min_precip_threshold=min_precip_threshold,
                                         logger=logger,
                                         cell_chunk_size=cell_chunk_size,
                                         available_memory_gb=available_memory_gb)
        
        # Compute results (convert from dask to numpy)
        logger.info("Computing results...")
        results_computed = {}
        for p, data in results.items():
            results_computed[p] = data.compute()
            logger.info(f"P{int(p)}: min={float(results_computed[p].min()):.3f}, "
                       f"max={float(results_computed[p].max()):.3f}, "
                       f"mean={float(results_computed[p].mean()):.3f} mm/h")
        
        # Create output filename
        time_str = time_duration.lower().replace('h', 'h').replace('d', 'd')
        out_basename = f"{source_name}_precip_percentiles_{time_str}_hp{zoom}_{version}.nc"
        output_filename = os.path.join(output_dir, out_basename)
        
        # Write to NetCDF
        write_netcdf(results_computed, ds_p, output_filename, zoom, source_name,
                    start_datetime, end_datetime, time_duration, method,
                    min_precip_threshold, logger, extra_attrs=extra_attrs)
    
    logger.info("=" * 60)
    logger.info("Extreme precipitation percentile calculation complete!")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
