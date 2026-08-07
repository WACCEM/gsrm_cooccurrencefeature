"""
Calculate extreme precipitation percentiles from HEALPix grid data (1-hourly).

This script computes precipitation percentiles (e.g., 90th, 95th) at different
time scales (1-hourly, 6-hourly, daily, etc.) from 1-hourly data on HEALPix grid.
All input data are sourced from the catalog, except IMERG which uses local 1-hourly
Zarr files.

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
import easygems.healpix as egh


def parse_cmd_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description='Calculate extreme precipitation percentiles (1-hourly input)')
    parser.add_argument('--catalog_source', type=str, required=True,
                        help='Catalog source name (e.g., scream_ne120, IR_IMERG)')
    parser.add_argument('--config_file', type=str,
                        default='/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources_1h.yaml',
                        help='Path to configuration YAML file '
                             '(default: /global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources_1h.yaml)')
    parser.add_argument('--zoom', type=int, default=8,
                        help='HEALPix zoom level (default: 8)')
    parser.add_argument('--percentiles', type=float, nargs='+', default=[95, 99],
                        help='Percentiles to compute (default: 95 99)')
    parser.add_argument('--time_durations', type=str, nargs='+', default=['1h'],
                        help='Time durations for resampling (default: 1h)')
    parser.add_argument('--method', type=str, default='linear',
                        choices=['linear', 'lower', 'higher', 'midpoint', 'nearest'],
                        help='Quantile interpolation method (default: linear)')
    parser.add_argument('--min_precip_threshold', type=float, default=0.1,
                        help='Minimum precipitation threshold (mm/h) to exclude from percentile calculation. '
                             'Values below this threshold are treated as missing. (default: 0.1 mm/h)')
    parser.add_argument('--min_year_coverage_days', type=int, default=300,
                        help='Minimum number of distinct calendar days with data required for a '
                             'calendar year to be included in the annual percentile calculation. '
                             '(default: 300)')
    parser.add_argument('--cell_chunk_size', type=int, default=None,
                        help='Number of HEALPix cells processed per block when computing '
                             'quantiles over time. Bounds peak memory by materializing one '
                             'block (all time steps x cell_chunk_size cells) at a time instead '
                             'of letting dask auto-rechunk the whole record into thousands of '
                             'tiny blocks, which is what causes out-of-memory failures on long '
                             '(multi-decade) hourly records. Default (None): sized '
                             'automatically per call -- as large as safely fits in available '
                             'memory, aligned to the native Zarr cell-chunk size to avoid '
                             'redundant chunk re-decompression (the main driver of runtime, not '
                             'just memory). Pass an explicit value to disable auto-sizing.')
    parser.add_argument('--available_memory_gb', type=float, default=None,
                        help='System memory (GB) to budget against when auto-sizing '
                             'cell_chunk_size. Default (None): auto-detected via psutil. Pass '
                             'this explicitly if running under a SLURM allocation where '
                             'node-level auto-detection may not reflect the job\'s actual '
                             'memory limit. Ignored if --cell_chunk_size is given explicitly.')
    parser.add_argument('--no_annual', action='store_true',
                        help='Skip computing annual percentiles and the interannual IQR '
                             '(default: annual/IQR are computed for sources with >= 2 '
                             'qualifying years)')
    parser.add_argument('--output_dir', type=str, default='/pscratch/sd/w/wcmca1/hackathon/extreme_precip_1h/',
                        help='Output directory for NetCDF files '
                             '(default: /pscratch/sd/w/wcmca1/hackathon/extreme_precip_1h/)')
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


def prepare_precip(pr, time_duration='1h', min_precip_threshold=None, logger=None):
    """
    Apply the minimum precipitation threshold and resample precipitation to the
    requested time duration. Shared preprocessing step used by both the all-period
    and the per-year percentile calculations, so it is only built once (lazily).

    Args:
        pr: xr.DataArray
            Precipitation data array with dimensions (time, cell) in mm/h
        time_duration: str
            Time duration for resampling. Options include:
            - '1h' or '1H': 1-hourly (no resampling if already 1-hourly)
            - '6h' or '6H': 6-hourly
            - '1D' or '24h': daily
            - '12h' or '12H': 12-hourly
            Default is '1h' (keeps original 1-hourly resolution)
        min_precip_threshold: float, optional
            Minimum precipitation threshold (mm/h). Values below this threshold
            are excluded from percentile calculations (treated as NaN).
            If None, all precipitation values are included. Default is None.
        logger: logging.Logger
            Logger instance

    Returns:
        xr.DataArray: Thresholded and resampled (still lazy) precipitation data
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Apply minimum precipitation threshold if specified
    if min_precip_threshold is not None:
        logger.info(f"Applying minimum precipitation threshold: {min_precip_threshold} mm/h")
        logger.info(f"Values below threshold will be excluded from percentile calculation")
        # Replace values below threshold with NaN
        pr_filtered = pr.where(pr >= min_precip_threshold)
    else:
        logger.info("No minimum precipitation threshold applied")
        pr_filtered = pr

    # Resample to specified time duration if needed.
    # Input data are 1-hourly; skip resampling when time_duration is also 1h.
    if time_duration.lower() not in ['1h', '1hr']:
        # Convert time duration to pandas-compatible frequency string
        freq = time_duration.upper().replace('HR', 'H')
        logger.info(f"Resampling precipitation from 1-hourly to {freq}...")
        # Resample by taking mean precipitation rate (ignoring NaN)
        pr_resampled = pr_filtered.resample(time=freq).mean(keep_attrs=True)
    else:
        logger.info(f"Using original 1-hourly precipitation data...")
        pr_resampled = pr_filtered

    return pr_resampled


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


def calc_precip_percentiles(pr_resampled, percentiles=[90, 95], time_duration='1h', method='linear',
                            min_precip_threshold=None, cell_chunk_size=None, available_memory_gb=None,
                            logger=None):
    """
    Calculate all-period precipitation percentiles from already-prepared data.

    Args:
        pr_resampled: xr.DataArray
            Precipitation data array with dimensions (time, cell) in mm/h, already
            passed through prepare_precip() (threshold applied, resampled to
            time_duration)
        percentiles: list of float
            Percentiles to compute (e.g., [90, 95] for 90th and 95th percentiles)
        time_duration: str
            Time duration used to prepare pr_resampled (only used here for
            variable attributes; see prepare_precip() for resampling options)
        method: str
            Interpolation method for quantile calculation. Options include:
            - 'linear': linear interpolation (default)
            - 'lower': lower value
            - 'higher': higher value
            - 'midpoint': midpoint of two nearest values
            - 'nearest': nearest value
        min_precip_threshold: float, optional
            Minimum precipitation threshold (mm/h) used to prepare pr_resampled
            (only used here for variable attributes). Default is None.
        cell_chunk_size: int, optional
            Number of cells processed per block; see determine_cell_chunk_size()
            for rationale. If None (default), sized dynamically.
        available_memory_gb: float, optional
            Forwarded to determine_cell_chunk_size() when cell_chunk_size is None.
        logger: logging.Logger
            Logger instance

    Returns:
        dict: Dictionary with percentile values as keys and DataArrays as values
              e.g., {90: pr_p90, 95: pr_p95}
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Calculate percentiles by looping over cell blocks to bound peak memory
    # (see _quantile_over_time_blocked for why this is necessary for multi-decade
    # hourly records)
    logger.info(f"Computing percentile(s) {percentiles} over the full record...")
    quantile_results = _quantile_over_time_blocked(pr_resampled, percentiles, method=method,
                                                    cell_chunk_size=cell_chunk_size,
                                                    available_memory_gb=available_memory_gb,
                                                    logger=logger)

    results = {}
    for p in percentiles:
        pr_percentile = quantile_results[p]
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


def calc_annual_precip_percentiles(pr_resampled, percentiles=[90, 95], time_duration='1h', method='linear',
                                   min_precip_threshold=None, min_year_coverage_days=300, min_years=2,
                                   cell_chunk_size=None, available_memory_gb=None, logger=None):
    """
    Calculate precipitation percentiles separately for each qualifying calendar year.

    Args:
        pr_resampled: xr.DataArray
            Precipitation data array with dimensions (time, cell) in mm/h, already
            passed through prepare_precip() (threshold applied, resampled to
            time_duration)
        percentiles: list of float
            Percentiles to compute (e.g., [90, 95] for 90th and 95th percentiles)
        time_duration: str
            Time duration used to prepare pr_resampled (only used here for
            variable attributes)
        method: str
            Quantile interpolation method (see calc_precip_percentiles)
        min_precip_threshold: float, optional
            Minimum precipitation threshold (mm/h) used to prepare pr_resampled
            (only used here for variable attributes). Default is None.
        min_year_coverage_days: int
            Minimum number of distinct calendar days with data required for a
            calendar year to be included in the annual calculation. Default is 300.
        min_years: int
            Minimum number of qualifying years required to compute annual
            percentiles at all. Default is 2.
        cell_chunk_size: int, optional
            Number of cells processed per block; see determine_cell_chunk_size()
            for rationale. If None (default), sized dynamically *per year* -- since
            a single year's time span is much shorter than the full record, this
            naturally lands on a much larger (often native-chunk-aligned, zero
            redundant-read) block size than the all-period call uses.
        available_memory_gb: float, optional
            Forwarded to determine_cell_chunk_size() when cell_chunk_size is None.
        logger: logging.Logger
            Logger instance

    Returns:
        tuple: (results, years)
            results: dict with percentile values as keys and DataArrays (year, cell)
                     as values, e.g., {90: pr_annual_p90, 95: pr_annual_p95}
            years: sorted list of qualifying calendar years (int)
            Returns (None, None) if fewer than min_years calendar years qualify.
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Determine data coverage (distinct calendar days present) for each calendar year.
    # Counting unique days (rather than time steps) keeps this independent of
    # time_duration and tolerant of gaps in the record.
    time_index = pd.DatetimeIndex(pr_resampled['time'].values)
    days = pd.Series(time_index.normalize()).drop_duplicates()
    coverage_days = days.groupby(days.dt.year).size()

    qualifying_years = sorted(int(y) for y, n in coverage_days.items() if n >= min_year_coverage_days)

    for y in coverage_days.index.sort_values():
        status = 'kept' if int(y) in qualifying_years else 'dropped (insufficient coverage)'
        logger.info(f"  Year {int(y)}: {coverage_days[y]} days of data -> {status}")

    if len(qualifying_years) < min_years:
        logger.info(f"Only {len(qualifying_years)} year(s) meet the minimum coverage of "
                    f"{min_year_coverage_days} days (need >= {min_years}); "
                    f"skipping annual percentile calculation.")
        return None, None

    logger.info(f"Computing annual percentiles for {len(qualifying_years)} year(s): {qualifying_years}")

    results = {p: [] for p in percentiles}
    for year in qualifying_years:
        logger.info(f"  Computing percentiles for year {year}...")
        pr_year = pr_resampled.sel(time=str(year))
        year_results = _quantile_over_time_blocked(pr_year, percentiles, method=method,
                                                    cell_chunk_size=cell_chunk_size,
                                                    available_memory_gb=available_memory_gb,
                                                    logger=logger)
        for p in percentiles:
            results[p].append(year_results[p])

    annual_results = {}
    for p in percentiles:
        pr_annual = xr.concat(results[p], dim=pd.Index(qualifying_years, name='year'))
        if 'quantile' in pr_annual.coords:
            pr_annual = pr_annual.drop_vars('quantile')
        pr_annual.name = f'pr_annual_p{p}'
        pr_annual.attrs['long_name'] = f'{p}th percentile precipitation for each calendar year'
        pr_annual.attrs['units'] = 'mm/h'
        pr_annual.attrs['time_duration'] = time_duration
        pr_annual.attrs['percentile'] = p
        pr_annual.attrs['method'] = method
        pr_annual.attrs['cell_methods'] = 'time: quantile (interval: 1 year)'
        pr_annual.attrs['comment'] = (f'Computed independently within each calendar year having '
                                       f'at least {min_year_coverage_days} days of data')
        pr_annual.attrs['min_year_coverage_days'] = min_year_coverage_days
        if min_precip_threshold is not None:
            pr_annual.attrs['min_precip_threshold'] = f'{min_precip_threshold} mm/h'
            pr_annual.attrs['note'] = f'Precipitation below {min_precip_threshold} mm/h excluded from calculation'
        annual_results[p] = pr_annual

    return annual_results, qualifying_years


def calc_interannual_iqr(annual_results, method='linear', logger=None):
    """
    Calculate the interquartile range (and quartiles) across years of the annual
    precipitation percentiles.

    Args:
        annual_results: dict
            Dictionary with percentile values as keys and computed DataArrays
            (dims: year, cell) as values, e.g., {90: pr_annual_p90, 95: pr_annual_p95}
        method: str
            Quantile interpolation method used for the q25/q75 calculation across years
        logger: logging.Logger
            Logger instance

    Returns:
        dict: Dictionary keyed by percentile, each value a dict with 'q25', 'q75',
              'iqr' DataArrays (dims: cell)
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    iqr_results = {}
    for p, pr_annual in annual_results.items():
        logger.info(f"Computing interannual IQR of the {p}th percentile across years...")
        n_years = pr_annual.sizes['year']
        quartiles = pr_annual.quantile([0.25, 0.75], dim='year', method=method, skipna=True)
        q25 = quartiles.sel(quantile=0.25, drop=True)
        q75 = quartiles.sel(quantile=0.75, drop=True)
        iqr = q75 - q25

        for name, da_, label in [('q25', q25, '25th'), ('q75', q75, '75th')]:
            da_.name = f'pr_{name}_p{p}'
            da_.attrs['long_name'] = f'{label} percentile across annual {p}th percentile precipitation'
            da_.attrs['units'] = 'mm/h'
            da_.attrs['percentile'] = p
            da_.attrs['method'] = method
            da_.attrs['cell_methods'] = 'year: quantile'
            da_.attrs['n_years'] = n_years
            da_.attrs['comment'] = (f'{label} percentile, taken across the {n_years} annual '
                                     f'{p}th percentile precipitation values, of the interannual '
                                     f'distribution of the {p}th percentile threshold')

        iqr.name = f'pr_iqr_p{p}'
        iqr.attrs['long_name'] = f'Interquartile range of annual {p}th percentile precipitation'
        iqr.attrs['units'] = 'mm/h'
        iqr.attrs['percentile'] = p
        iqr.attrs['method'] = method
        iqr.attrs['cell_methods'] = 'year: quantile'
        iqr.attrs['n_years'] = n_years
        iqr.attrs['comment'] = (f'pr_q75_p{p} minus pr_q25_p{p}; measures the interannual spread '
                                 f'of the {p}th percentile precipitation threshold across {n_years} years')

        iqr_results[p] = {'q25': q25, 'q75': q75, 'iqr': iqr}

    return iqr_results


def calc_snow_probability_mean(ds_p, varname='snowProbability', logger=None):
    """
    Calculate the time-mean snow probability, if the variable is present in the dataset.

    This provides backward compatibility for sources that do not carry a snow
    probability field (e.g. IMERG, and the catalog model sources): the function simply
    returns None, and callers skip adding the variable to the output.

    Args:
        ds_p: xr.Dataset
            Source dataset (as returned by load_precipitation_data), potentially
            containing a snow probability variable with dimensions (time, cell)
        varname: str
            Name of the snow probability variable to look for. Default is
            'snowProbability' (as found in the GSMaP HEALPix Zarr files).
        logger: logging.Logger
            Logger instance

    Returns:
        xr.DataArray or None: Time-mean snow probability (dim: cell), still lazy, or
            None if `varname` is not present in ds_p (or lacks a time dimension)
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    if varname not in ds_p:
        logger.info(f"Snow probability variable '{varname}' not found in dataset; skipping.")
        return None

    snow_prob_in = ds_p[varname]
    if 'time' not in snow_prob_in.dims:
        logger.info(f"'{varname}' found but has no time dimension; skipping.")
        return None

    logger.info(f"Found snow probability variable '{varname}'; computing time mean...")
    units = snow_prob_in.attrs.get('units', 'percent')

    snow_prob = snow_prob_in.mean(dim='time', skipna=True)
    snow_prob.name = 'snow_probability'
    snow_prob.attrs['long_name'] = 'Time-mean probability of snow'
    snow_prob.attrs['units'] = units
    snow_prob.attrs['cell_methods'] = 'time: mean'
    snow_prob.attrs['source_variable'] = varname
    snow_prob.attrs['comment'] = (f'Mean of {varname} over all times in the record; '
                                   f'not filtered by min_precip_threshold and independent '
                                   f'of the requested time_duration')

    return snow_prob


def write_netcdf(results_dict, ds_p, output_filename, zoom, source_name,
                 start_datetime, end_datetime, time_duration, method,
                 min_precip_threshold=None, logger=None,
                 annual_results=None, iqr_results=None, years=None,
                 min_year_coverage_days=None, snow_probability=None):
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
        annual_results: dict, optional
            Dictionary with percentile values as keys and computed DataArrays
            (dims: year, cell) as values. If None, no annual variables are written.
        iqr_results: dict, optional
            Dictionary keyed by percentile, each value a dict with 'q25', 'q75',
            'iqr' computed DataArrays (dims: cell). If None, no interannual IQR
            variables are written.
        years: list of int, optional
            Qualifying calendar years corresponding to annual_results' year dimension
        min_year_coverage_days: int, optional
            Minimum coverage (days) required for a year to qualify, recorded as an
            attribute
        snow_probability: xr.DataArray, optional
            Computed time-mean snow probability (dim: cell), as returned by
            calc_snow_probability_mean(). If None, no snow_probability variable is
            written (backward compatible for sources without the field).
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

    # Add annual percentile and interannual IQR variables, if provided
    if annual_results is not None and years is not None:
        coord_dict['year'] = (['year'], np.asarray(years, dtype='int32'))
        for p, data in annual_results.items():
            var_name = f'pr_annual_p{int(p)}'
            var_dict[var_name] = (['year', 'cell'], data.values)

    if iqr_results is not None:
        for p, quartiles in iqr_results.items():
            for name in ('q25', 'q75', 'iqr'):
                var_name = f'pr_{name}_p{int(p)}'
                var_dict[var_name] = (['cell'], quartiles[name].values)

    # Add time-mean snow probability, if computed (absent for sources without the field)
    if snow_probability is not None:
        var_dict['snow_probability'] = (['cell'], snow_probability.values)

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

    if years is not None:
        gattr_dict['annual_years'] = ', '.join(str(y) for y in years)
        gattr_dict['n_annual_years'] = len(years)
        if min_year_coverage_days is not None:
            gattr_dict['min_year_coverage_days'] = min_year_coverage_days

    if snow_probability is not None:
        gattr_dict['snow_probability_source'] = snow_probability.attrs.get('source_variable', '')

    # Create output dataset
    dsout = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)

    # Add coordinate attributes
    dsout['cell'].attrs['long_name'] = 'HEALPix cell index'
    dsout['lon'].attrs['long_name'] = 'Longitude'
    dsout['lon'].attrs['units'] = 'degree'
    dsout['lat'].attrs['long_name'] = 'Latitude'
    dsout['lat'].attrs['units'] = 'degree'
    if 'year' in dsout.coords:
        dsout['year'].attrs['long_name'] = 'Calendar year'
        dsout['year'].attrs['description'] = (
            f'Calendar years with at least {min_year_coverage_days} days of data, '
            f'used for the annual percentile and interannual IQR variables'
            if min_year_coverage_days is not None else
            'Calendar years used for the annual percentile and interannual IQR variables'
        )

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

    # Carry over attributes computed on the source DataArrays for the annual and
    # interannual IQR variables (long_name, units, cell_methods, comment, etc.)
    if annual_results is not None and years is not None:
        for p, data in annual_results.items():
            var_name = f'pr_annual_p{int(p)}'
            dsout[var_name].attrs.update(data.attrs)

    if iqr_results is not None:
        for p, quartiles in iqr_results.items():
            for name in ('q25', 'q75', 'iqr'):
                var_name = f'pr_{name}_p{int(p)}'
                dsout[var_name].attrs.update(quartiles[name].attrs)
                if min_precip_threshold is not None:
                    dsout[var_name].attrs['min_precip_threshold'] = f'{min_precip_threshold} mm/h'

    if snow_probability is not None:
        dsout['snow_probability'].attrs.update(snow_probability.attrs)

    # Save the output file
    fillvalue = np.nan
    comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}

    logger.info(f'Writing output file: {output_filename}')
    dsout.to_netcdf(path=output_filename, mode='w', format='NETCDF4', encoding=encoding)
    logger.info(f'Successfully wrote: {output_filename}')

    return dsout


def load_precipitation_data(config_file, catalog_source, zoom, logger=None):
    """
    Load 1-hourly precipitation data from catalog or direct Zarr file.

    All data sources are loaded from the HEALPix catalog, with the exception
    of IMERG which uses local 1-hourly Zarr files.

    Args:
        config_file: str
            Path to configuration YAML file
        catalog_source: str
            Catalog source name
        zoom: int
            HEALPix zoom level
        logger: logging.Logger
            Logger instance

    Returns:
        tuple: (pr DataArray, ds_p Dataset, config dict)
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

    # Catalog file
    catalog_file = "https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml"

    if catalog_source == "IMERG":
        # IMERG 1-hourly data is not in the catalog; use local Zarr files
        dir_healpix = "/pscratch/sd/w/wcmca1/GPM/healpix/"
        # in_zarr = f"{dir_healpix}IMERG_V7_1H_zoom{zoom}_20190101_20211231.zarr"
        in_zarr = f"{dir_healpix}IMERG_V7_1H_zoom{zoom}_20010101_20241231.zarr"
        logger.info(f"Loading IMERG 1-hourly dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)

    elif catalog_source == "GSMAP":
        # GSMAP 1-hourly data is not in the catalog; use local Zarr files
        dir_healpix = "/pscratch/sd/w/wcmca1/GsMAP/healpix/"
        in_zarr = f"{dir_healpix}GsMAPv8_1H_zoom{zoom}_20100101_20241231.zarr"
        # in_zarr = f"{dir_healpix}GsMAPv8_1H_zoom{zoom}_20200101_20201231.zarr"
        logger.info(f"Loading GSMAP 1-hourly dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)

    else:
        # All other sources: load from the HEALPix catalog
        # Update the zoom level in catalog_params
        catalog_params['zoom'] = zoom
        logger.info(f"Loading HEALPix catalog: {catalog_file}")
        in_catalog = intake.open_catalog(catalog_file)
        if catalog_location:
            in_catalog = in_catalog[catalog_location]

        logger.info(f"Loading {catalog_source} from catalog with params: {catalog_params}")
        ds_p = in_catalog[catalog_source](**catalog_params).to_dask()
        # Add lat/lon coordinates to the HEALPix DataSet
        ds_p = ds_p.pipe(egh.attach_coords)
    # import pdb; pdb.set_trace()

    # Check liquid precipitation variable
    if varname_precip_liq in list(ds_p.keys()):
        # Convert liquid precipitation to mm/h
        pr = ds_p[varname_precip_liq] * pr_convert_factor

    # Check if the ice precipitation variable exists in the dataset
    if varname_precip_ice and varname_precip_ice in list(ds_p.keys()):
        # Convert ice precipitation to liquid equivalent
        prs = ds_p[varname_precip_ice] * pr_convert_factor
        # Add ice precipitation to get total precipitation
        pr = pr + prs

    logger.info(f"Precipitation data loaded: {pr.shape}")

    return pr, ds_p, config


def main():
    """
    Main function to calculate extreme precipitation percentiles from 1-hourly data.
    """
    # Parse command-line arguments
    args_dict = parse_cmd_args()

    # Setup logging
    logger = setup_logging()
    logger.info("=" * 60)
    logger.info("Starting extreme precipitation percentile calculation (1-hourly)")
    logger.info("=" * 60)

    # Extract arguments
    config_file = args_dict['config_file']
    catalog_source = args_dict['catalog_source']
    zoom = args_dict['zoom']
    percentiles = args_dict['percentiles']
    time_durations = args_dict['time_durations']
    method = args_dict['method']
    min_precip_threshold = args_dict['min_precip_threshold']
    min_year_coverage_days = args_dict['min_year_coverage_days']
    cell_chunk_size = args_dict['cell_chunk_size']
    available_memory_gb = args_dict['available_memory_gb']
    compute_annual = not args_dict['no_annual']
    output_dir = args_dict['output_dir']
    version = args_dict['version']

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Load precipitation data
    logger.info(f"Loading 1-hourly precipitation data for {catalog_source}...")
    pr, ds_p, config = load_precipitation_data(config_file, catalog_source, zoom, logger)

    source_name = config.get('source_name')
    start_datetime = config.get('start_datetime')
    end_datetime = config.get('end_datetime')
    varname_snow_probability = config.get('varname_snow_probability', 'snowProbability')

    # Calculate time-mean snow probability once (independent of time_duration).
    # Returns None (and is skipped everywhere downstream) if the source has no
    # snow probability variable, keeping this backward compatible.
    snow_prob_computed = calc_snow_probability_mean(ds_p, varname=varname_snow_probability, logger=logger)
    if snow_prob_computed is not None:
        logger.info("Computing time-mean snow probability...")
        snow_prob_computed = snow_prob_computed.compute()
        logger.info(f"Snow probability: min={float(snow_prob_computed.min()):.3f}, "
                    f"max={float(snow_prob_computed.max()):.3f}, "
                    f"mean={float(snow_prob_computed.mean()):.3f} {snow_prob_computed.attrs.get('units', '')}")

    # Process each time duration
    for time_duration in time_durations:
        logger.info("=" * 60)
        logger.info(f"Processing time duration: {time_duration}")
        logger.info("=" * 60)

        # Apply threshold and resample once; reused by both the all-period and
        # annual percentile calculations
        pr_resampled = prepare_precip(pr, time_duration=time_duration,
                                      min_precip_threshold=min_precip_threshold,
                                      logger=logger)

        # Calculate all-period percentiles
        results = calc_precip_percentiles(pr_resampled, percentiles=percentiles,
                                          time_duration=time_duration,
                                          method=method,
                                          min_precip_threshold=min_precip_threshold,
                                          cell_chunk_size=cell_chunk_size,
                                          available_memory_gb=available_memory_gb,
                                          logger=logger)

        # Compute results (convert from dask to numpy)
        logger.info("Computing results...")
        results_computed = {}
        for p, data in results.items():
            results_computed[p] = data.compute()
            logger.info(f"P{int(p)}: min={float(results_computed[p].min()):.3f}, "
                        f"max={float(results_computed[p].max()):.3f}, "
                        f"mean={float(results_computed[p].mean()):.3f} mm/h")

        # Calculate annual percentiles and their interannual IQR, for sources with
        # enough qualifying calendar years
        annual_computed = None
        iqr_computed = None
        years = None
        if compute_annual:
            logger.info("-" * 60)
            logger.info("Checking calendar year coverage for annual percentiles...")
            annual_results, years = calc_annual_precip_percentiles(
                pr_resampled, percentiles=percentiles, time_duration=time_duration,
                method=method, min_precip_threshold=min_precip_threshold,
                min_year_coverage_days=min_year_coverage_days,
                cell_chunk_size=cell_chunk_size, available_memory_gb=available_memory_gb,
                logger=logger)

            if annual_results is not None:
                logger.info("Computing annual results...")
                annual_computed = {}
                for p, data in annual_results.items():
                    annual_computed[p] = data.compute()
                    for iy, year in enumerate(years):
                        vals = annual_computed[p].isel(year=iy)
                        logger.info(f"P{int(p)} {year}: min={float(vals.min()):.3f}, "
                                    f"max={float(vals.max()):.3f}, "
                                    f"mean={float(vals.mean()):.3f} mm/h")

                logger.info("Computing interannual IQR...")
                iqr_computed = calc_interannual_iqr(annual_computed, method=method, logger=logger)
                for p, quartiles in iqr_computed.items():
                    iqr = quartiles['iqr']
                    logger.info(f"IQR P{int(p)}: min={float(iqr.min()):.3f}, "
                                f"max={float(iqr.max()):.3f}, "
                                f"mean={float(iqr.mean()):.3f} mm/h")

        # Create output filename
        time_str = time_duration.lower().replace('h', 'h').replace('d', 'd')
        out_basename = f"{source_name}_precip_percentiles_{time_str}_hp{zoom}_{version}.nc"
        output_filename = os.path.join(output_dir, out_basename)

        # Write to NetCDF
        write_netcdf(results_computed, ds_p, output_filename, zoom, source_name,
                     start_datetime, end_datetime, time_duration, method,
                     min_precip_threshold, logger,
                     annual_results=annual_computed, iqr_results=iqr_computed,
                     years=years, min_year_coverage_days=min_year_coverage_days,
                     snow_probability=snow_prob_computed)

    logger.info("=" * 60)
    logger.info("Extreme precipitation percentile calculation complete!")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
