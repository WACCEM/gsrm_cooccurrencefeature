import xarray as xr
import numpy as np
import pandas as pd
import cftime
import yaml
import calendar
import os, glob, re
import time
import argparse
import logging
import traceback
import sys
import gc
import intake
import easygems.healpix as egh
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from src.zarr_tools import setup_dask_client, initialize_zarr_store, append_chunk_to_zarr
from src.utilities import convert_cftime_to_standard_calendar
from pyflextrkr.ft_utilities import load_config
from pyflextrkr.ftfunctions import olr_to_tb

# Import for parallel processing
try:
    from distributed import as_completed
except ImportError:
    as_completed = None  # Will only be needed if parallel=True

def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Suppress verbose Dask logging
    logging.getLogger('distributed').setLevel(logging.WARNING)
    logging.getLogger('distributed.worker').setLevel(logging.WARNING)
    logging.getLogger('distributed.core').setLevel(logging.WARNING)
    logging.getLogger('distributed.comm').setLevel(logging.WARNING)
    logging.getLogger('distributed.nanny').setLevel(logging.WARNING)
    logging.getLogger('distributed.scheduler').setLevel(logging.WARNING)

#--------------------------------------------------------------------------------------------------
def create_track_swaths_and_coverage(tracknumber):
    """
    Create track swaths and coverage count arrays from a 3D track number array.
    
    Parameters:
    -----------
    tracknumber : numpy.ndarray
        3D array with dimensions (time, y, x) containing track numbers at each timestep.
        Background/no-track pixels should be 0.
    
    Returns:
    --------
    track_swaths_dict : dict
        Dictionary with track_id as keys and 2D swath arrays as values.
        Each swath shows all pixels covered by that track (labeled with track_id).
    track_coverage_dict : dict
        Dictionary with track_id as keys and 2D coverage count arrays as values.
        Each array shows the number of timesteps each pixel was covered by that track.
    
    Example:
    --------
    >>> tracknumber = np.zeros((24, 200, 150), dtype=int)
    >>> # ... populate tracknumber with track data ...
    >>> swaths, coverage = create_track_swaths_and_coverage(tracknumber)
    """
    # Get unique track numbers (excluding 0 which is background)
    unique_tracks = np.unique(tracknumber)
    unique_tracks = unique_tracks[unique_tracks > 0]
    
    # Initialize dictionaries
    track_swaths_dict = {}
    track_coverage_dict = {}
    
    # Process each unique track using vectorized operations
    for track_id in unique_tracks:
        # Create binary mask for this track across all timesteps
        track_mask = (tracknumber == track_id).astype(int)
        
        # Track swath: all pixels covered by this track labeled with track number
        swath = np.any(track_mask, axis=0).astype(int) * track_id
        track_swaths_dict[track_id] = swath
        
        # Coverage count: number of times each pixel is covered
        coverage = np.sum(track_mask, axis=0)
        track_coverage_dict[track_id] = coverage
    
    return track_swaths_dict, track_coverage_dict

#--------------------------------------------------------------------------------------------------
def combine_swaths_with_priority(track_swaths_dict, track_coverage_dict):
    """
    Combine multiple track swaths into a single 2D array, resolving overlaps
    by assigning the track number with the largest coverage count.
    
    Parameters:
    -----------
    track_swaths_dict : dict
        Dictionary with track_id as keys and 2D swath arrays as values
    track_coverage_dict : dict
        Dictionary with track_id as keys and 2D coverage count arrays as values
    
    Returns:
    --------
    combined_swath : numpy.ndarray
        2D array with combined swaths, where overlapping pixels are assigned
        to the track with the highest coverage count
    """
    # Get the shape from the first swath
    first_key = list(track_swaths_dict.keys())[0]
    shape = track_swaths_dict[first_key].shape
    
    # Initialize output array
    combined_swath = np.zeros(shape, dtype=int)
    
    # Initialize array to track maximum coverage at each pixel
    max_coverage = np.zeros(shape, dtype=int)
    
    # Process each track using vectorized operations
    for track_id in track_swaths_dict.keys():
        # Get swath and coverage for this track
        swath = track_swaths_dict[track_id]
        coverage = track_coverage_dict[track_id]
        
        # Create mask for pixels belonging to this track's swath
        track_mask = swath > 0
        
        # Update combined_swath where:
        # 1. Current pixel has no assignment yet (combined_swath == 0), OR
        # 2. Current track has higher coverage than previous assignment
        update_mask = track_mask & ((combined_swath == 0) | (coverage > max_coverage))
        
        # Apply updates using vectorized operations
        combined_swath = np.where(update_mask, track_id, combined_swath)
        max_coverage = np.where(update_mask, coverage, max_coverage)
    
    return combined_swath

#--------------------------------------------------------------------------------------------------
def create_latitude_dependent_tb_threshold(lat_values):
    """
    Create latitude-dependent brightness temperature thresholds.
    
    Parameters:
    -----------
    lat_values : numpy.ndarray
        Array of latitude values
        
    Returns:
    --------
    tb_thresh : numpy.ndarray
        Array of temperature thresholds (K) matching lat_values shape
    """
    abs_lat = np.abs(lat_values)
    tb_thresh = np.zeros_like(lat_values, dtype=np.float32)
    
    # Tropics: |lat| <= 30, tb_thresh = 250 K
    tropical_mask = abs_lat <= 30
    tb_thresh[tropical_mask] = 250.0
    
    # Mid-latitudes: 30 < |lat| <= 60, linearly decrease from 250 to 230 K
    midlat_mask = (abs_lat > 30) & (abs_lat <= 60)
    tb_thresh[midlat_mask] = 250.0 - 20.0 * ((abs_lat[midlat_mask] - 30.0) / 30.0)
    
    # High latitudes: |lat| > 60, tb_thresh = 230 K
    highlat_mask = abs_lat > 60
    tb_thresh[highlat_mask] = 230.0
    
    return tb_thresh

#--------------------------------------------------------------------------------------------------
def classify_cloud_types(tb, pr, tb_thresh, mcs_mask, pr_threshold=0.5):
    """
    Classify cloud types based on brightness temperature and precipitation.
    
    Parameters:
    -----------
    tb : numpy.ndarray
        Brightness temperature (K)
    pr : numpy.ndarray
        Precipitation rate (mm/h)
    tb_thresh : numpy.ndarray
        Latitude-dependent brightness temperature threshold (K)
    mcs_mask : numpy.ndarray
        MCS mask array (values > 0 indicate MCS pixels)
        Classification only occurs where mcs_mask == 0 (non-MCS areas)
    pr_threshold : float, optional
        Precipitation threshold for classification (default: 0.5 mm/h)
    
    Returns:
    --------
    cloud_type : numpy.ndarray
        Cloud type classification:
        0 = Unclassified (inside MCS regions)
        1 = Deep convective (tb < tb_thresh & pr >= pr_threshold)
        2 = Stratiform (tb < tb_thresh & pr < pr_threshold)
        3 = Non-deep convective (tb >= tb_thresh & pr >= pr_threshold)
        4 = Drizzle (tb >= tb_thresh & pr < pr_threshold)
    """
    # Initialize cloud type array with zeros
    cloud_type = np.zeros_like(tb, dtype=np.int8)
    
    # Only classify where mcs_mask == 0 (non-MCS areas)
    valid_mask = mcs_mask == 0
    
    # Classification conditions (only applied where valid_mask is True)
    # 1. Deep convective: tb < tb_thresh & pr >= pr_threshold
    deep_conv = (tb < tb_thresh) & (pr >= pr_threshold) & valid_mask
    cloud_type[deep_conv] = 1
    
    # 2. Stratiform: tb < tb_thresh & pr < pr_threshold
    stratiform = (tb < tb_thresh) & (pr < pr_threshold) & valid_mask
    cloud_type[stratiform] = 2
    
    # 3. Non-deep convective: tb >= tb_thresh & pr >= pr_threshold
    nondeep_conv = (tb >= tb_thresh) & (pr >= pr_threshold) & valid_mask
    cloud_type[nondeep_conv] = 3
    
    # 4. Drizzle: tb >= tb_thresh & pr < pr_threshold
    drizzle = (tb >= tb_thresh) & (pr < pr_threshold) & valid_mask
    cloud_type[drizzle] = 4
    
    return cloud_type

#--------------------------------------------------------------------------------------------------
def find_priority_based_cloud_type(cloud_types_time_series):
    """
    Find cloud type based on priority ranking: 1 > 2 > 3 > 4 (VECTORIZED).
    
    If a cell has any occurrence of type 1, assign type 1.
    If no type 1 but has type 2, assign type 2.
    If no type 1 or 2 but has type 3, assign type 3.
    If no type 1, 2, or 3 but has type 4, assign type 4.
    Otherwise, assign 0 (unclassified).
    
    Parameters:
    -----------
    cloud_types_time_series : numpy.ndarray
        Array of shape (n_times, n_cells) with cloud type values starting from 0
        
    Returns:
    --------
    priority_based : numpy.ndarray
        Array of shape (n_cells,) with priority-based cloud type
    """
    n_times, n_cells = cloud_types_time_series.shape
    
    # Initialize result array with zeros
    priority_based = np.zeros(n_cells, dtype=np.int8)
    
    # Check for presence of each cloud type (any occurrence over time)
    has_type_1 = np.any(cloud_types_time_series == 1, axis=0)
    has_type_2 = np.any(cloud_types_time_series == 2, axis=0)
    has_type_3 = np.any(cloud_types_time_series == 3, axis=0)
    has_type_4 = np.any(cloud_types_time_series == 4, axis=0)
    
    # Apply priority: 1 > 2 > 3 > 4
    # Start from lowest priority and work up (so higher priorities overwrite)
    priority_based[has_type_4] = 4
    priority_based[has_type_3] = 3
    priority_based[has_type_2] = 2
    priority_based[has_type_1] = 1
    
    return priority_based

#--------------------------------------------------------------------------------------------------
def compute_mean_precip_by_cloud_type(cloud_types_timeseries, pr_timeseries):
    """
    Compute frequency-weighted mean precipitation for each cloud type over the aggregation window.
    
    This computes: (conditional mean) × (frequency of that type)
    When summed across all 4 types, this equals the simple mean of all precipitation.

    Parameters
    ----------
    cloud_types_timeseries : np.ndarray
        Array of shape (n_times, n_cells) with cloud type values (0-4).
    pr_timeseries : np.ndarray
        Array of shape (n_times, n_cells) with precipitation values (mm/h).

    Returns
    -------
    dc_pr : np.ndarray
        Frequency-weighted deep convective precipitation (cloud_type==1), shape (n_cells,)
    st_pr : np.ndarray
        Frequency-weighted stratiform precipitation (cloud_type==2), shape (n_cells,)
    nd_pr : np.ndarray
        Frequency-weighted non-deep convective precipitation (cloud_type==3), shape (n_cells,)
    dz_pr : np.ndarray
        Frequency-weighted drizzle precipitation (cloud_type==4), shape (n_cells,)
        
    Note
    ----
    dc_pr + st_pr + nd_pr + dz_pr = mean(pr_timeseries, axis=0)
    """
    import warnings
    
    n_times = cloud_types_timeseries.shape[0]
    
    # Mask precipitation by cloud type for each time
    dc_mask = cloud_types_timeseries == 1
    st_mask = cloud_types_timeseries == 2
    nd_mask = cloud_types_timeseries == 3
    dz_mask = cloud_types_timeseries == 4

    # Compute conditional means (average precipitation WHEN cell is this type)
    # Suppress RuntimeWarnings for empty slices
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', r'Mean of empty slice')
        dc_pr_cond = np.nanmean(np.where(dc_mask, pr_timeseries, np.nan), axis=0)
        st_pr_cond = np.nanmean(np.where(st_mask, pr_timeseries, np.nan), axis=0)
        nd_pr_cond = np.nanmean(np.where(nd_mask, pr_timeseries, np.nan), axis=0)
        dz_pr_cond = np.nanmean(np.where(dz_mask, pr_timeseries, np.nan), axis=0)
    
    # Compute frequency of each type (fraction of time steps)
    dc_freq = dc_mask.sum(axis=0) / n_times
    st_freq = st_mask.sum(axis=0) / n_times
    nd_freq = nd_mask.sum(axis=0) / n_times
    dz_freq = dz_mask.sum(axis=0) / n_times
    
    # Weight conditional means by frequency
    # Set NaN conditional means to 0 before weighting (cells that never had this type)
    dc_pr = np.nan_to_num(dc_pr_cond, nan=0.0) * dc_freq
    st_pr = np.nan_to_num(st_pr_cond, nan=0.0) * st_freq
    nd_pr = np.nan_to_num(nd_pr_cond, nan=0.0) * nd_freq
    dz_pr = np.nan_to_num(dz_pr_cond, nan=0.0) * dz_freq

    return dc_pr, st_pr, nd_pr, dz_pr

#--------------------------------------------------------------------------------------------------
def add_tb_pr_to_dataset(_ds, config):
    """
    Add brightness temperature (tb) and precipitation (pr) variables to a dataset chunk.
    
    This function reads OLR and precipitation data from a catalog, handles calendar conversions
    if necessary, and adds tb and pr variables to the input dataset.
    
    Parameters:
    -----------
    _ds : xarray.Dataset
        Input dataset chunk with time coordinate
    config : dict
        Configuration dictionary containing catalog information and variable names:
        - catalog_file: Path to the intake catalog file
        - catalog_location: Location within the catalog
        - catalog_source: Source name within the catalog
        - catalog_params: Parameters for the catalog source
        - varname_precip_liq: Variable name for liquid precipitation
        - varname_precip_ice: Variable name for ice precipitation
        - pcp_convert_factor: Conversion factor for precipitation (e.g., to convert to mm/h)
    
    Returns:
    --------
    _ds : xarray.Dataset
        Dataset with added 'tb' and 'pr' variables
    """
    # Extract config values
    catalog_file = config['catalog_file']
    catalog_location = config['catalog_location']
    catalog_source = config['catalog_source']
    catalog_params = config['catalog_params']
    varname_olr = config['varname_olr']
    varname_precip_liq = config['varname_precip_liq']
    varname_precip_ice = config['varname_precip_ice']
    pr_convert_factor = config['pcp_convert_factor']
    
    # Special handling for IR_IMERG: tb is already in the dataset, no OLR conversion needed
    if catalog_source == 'IR_IMERG':
        varname_olr = 'Tb'
        varname_precip_liq = 'precipitation'

    # Read OLR/precipitation data from catalog
    if catalog_location is None:
        cat = intake.open_catalog(catalog_file)
    else:
        cat = intake.open_catalog(catalog_file)[catalog_location]
    ds_p = cat[catalog_source](**catalog_params).to_dask()

    # List of variables to keep
    vars_keep = [varname_olr, varname_precip_liq, varname_precip_ice]
    # Identify variables to drop (those not in the keep list)
    dropvars_list = [v for v in ds_p.data_vars if v not in vars_keep]
    # Drop unwanted variables
    ds_p = ds_p.drop_vars(dropvars_list, errors='ignore')

    # 1. Determine calendar types and convert if needed
    ds_p_calendar_type = type(ds_p.time.values[0]).__name__
    ds_calendar_type = type(_ds.time.values[0]).__name__

    # Convert ds_p time to match _ds if they differ (convert non-standard to standard)
    if ds_p_calendar_type != ds_calendar_type:
        # Convert cftime objects to numpy datetime64 (standard calendar)
        converted_times = convert_cftime_to_standard_calendar(ds_p.time.values)
        
        # Create new dataset with converted time coordinate
        ds_p = ds_p.assign_coords(time=converted_times)

    # 2. Find common time range between datasets
    common_times = sorted(set(ds_p['time'].values)
                            .intersection(set(_ds['time'].values)))
    
    if not common_times:
        raise ValueError("No common time values between mask dataset and catalog dataset!")
    
    # 3. Select ds_p to only the common times (BEFORE expensive operations).
    # Eagerly load the small time slice into memory here so the full multi-year
    # catalog task graph (with its native 262144-cell chunks) is released before
    # we derive pr/tb. Without .load(), the 4x chunk size mismatch between the
    # catalog (cell=262144) and the MCS mask zarr (cell=65536) causes a large
    # intermediate memory spike when _ds.load() is called later in the worker.
    ds_p = ds_p.sel(time=common_times).load()
    
    # 4. Convert precipitation units (only for common times)
    # Check liquid precipitation variable
    if varname_precip_liq in list(ds_p.keys()):
        # Convert liquid precipitation to mm/h (if conversion factor is provided)
        if pr_convert_factor is not None:
            pr = ds_p[varname_precip_liq] * pr_convert_factor
        else:
            pr = ds_p[varname_precip_liq]
    # Check if the ice precipitation variable exist in the dataset
    if varname_precip_ice in list(ds_p.keys()):
        # Convert ice precipitation to liquid equivalent (if conversion factor is provided)
        if pr_convert_factor is not None:
            prs = ds_p[varname_precip_ice] * pr_convert_factor
        else:
            prs = ds_p[varname_precip_ice]
        # Add ice precipitation to get total precipitation
        pr = pr + prs
    
    # Convert OLR to Tb (only for common times)
    # For IR_IMERG, tb is already in the dataset, skip conversion
    if catalog_source == 'IR_IMERG':
        tb = ds_p[varname_olr]  # varname_olr is already set to 'tb' for IR_IMERG
    else:
        tb = olr_to_tb(ds_p[varname_olr])
    
    # 5. Add precipitation & tb to the dataset
    _ds = _ds.sel(time=common_times)
    _ds["pr"] = pr
    _ds["tb"] = tb
    
    return _ds

#--------------------------------------------------------------------------------------------------
def process_timechunk_swath(_ds, tb_thresh=None, verbose=False):
    """
    Process a time chunk of dataset to create MCS swath masks and cloud type classification.
    
    Parameters:
    -----------
    _ds : xarray.Dataset
        Input dataset chunk with dimensions (time, cell) containing 'mcs_mask', 'ccs_mask', 
        'tb', 'pr', and coordinates 'lat', 'lon'.
    tb_thresh : numpy.ndarray, optional
        Pre-computed latitude-dependent brightness temperature threshold.
        If None, will be computed from lat coordinate.
    verbose : bool
        If True, print progress information.
    
    Returns:
    --------
    dict : Dictionary containing:
        'mcs_mask': MCS swath mask (1D array, cell dimension)
        'ccs_mask': Non-MCS CCS mask (1D array, cell dimension)
        'cloud_types': Aggregated cloud type classification (1D array, cell dimension)
    """
    if verbose:
        print(f"Processing time chunk with {len(_ds.time)} time steps...")
    
    # Extract track number arrays
    mcs_mask = _ds['mcs_mask'].values  # shape (time, cell)
    # Replace NaN with 0 (when mask_and_scale=True, fill_value becomes NaN)
    mcs_mask = np.nan_to_num(mcs_mask, nan=0.0).astype(int)
    
    # # Sum CCS mask over time and convert to binary
    # # First replace NaN with 0, then sum
    # ccs_mask_values = _ds['ccs_mask'].values
    # ccs_mask_values = np.nan_to_num(ccs_mask_values, nan=0.0)
    # ccs_mask_sum = ((ccs_mask_values > 0).sum(axis=0) > 0).astype(int)  # shape (cell)
    
    # Create swaths and coverage for MCS
    mcs_swaths_dict, mcs_coverage_dict = create_track_swaths_and_coverage(mcs_mask)
    combined_mcs_swath = combine_swaths_with_priority(mcs_swaths_dict, mcs_coverage_dict)

    # Filter out CCS that overlap with MCS swaths
    # ccs_mask_sum[combined_mcs_swath > 0] = 0
    
    # ===== Cloud Type Classification =====
    # Create latitude-dependent tb threshold if not provided
    if tb_thresh is None:
        lat_values = _ds['lat'].values
        tb_thresh = create_latitude_dependent_tb_threshold(lat_values)
    
    # Classify cloud types for each time step
    cloud_types_timeseries = []
    for t in range(len(_ds.time)):
        tb_t = _ds['tb'].isel(time=t).values
        pr_t = _ds['pr'].isel(time=t).values
        mcs_mask_t = _ds['mcs_mask'].isel(time=t).values
        
        # Replace NaN with 0 (when mask_and_scale=True, fill_value becomes NaN)
        mcs_mask_t = np.nan_to_num(mcs_mask_t, nan=0.0).astype(int)
        
        # Classify cloud types (only in non-MCS areas)
        cloud_type_t = classify_cloud_types(tb_t, pr_t, tb_thresh, mcs_mask_t, pr_threshold=0.5)
        cloud_types_timeseries.append(cloud_type_t)
    
    # Stack into array (time, cell)
    cloud_types_timeseries = np.stack(cloud_types_timeseries, axis=0)
    
    # Extract precipitation time series for computing precipitation by cloud type
    pr_timeseries = _ds['pr'].values  # shape (time, cell)
    
    # Compute frequency-weighted mean precipitation for each cloud type
    dc_pr, st_pr, nd_pr, dz_pr = compute_mean_precip_by_cloud_type(cloud_types_timeseries, pr_timeseries)
    
    # Find cloud type based on priority ranking (1 > 2 > 3 > 4)
    cloud_types_aggregated = find_priority_based_cloud_type(cloud_types_timeseries)
    
    if verbose:
        n_classified_before = (cloud_types_aggregated > 0).sum()
        print(f"  Cloud types before MCS filter: {n_classified_before:,} cells")
    
    # Apply MCS priority: Set cloud type to 0 where MCS swath exists
    # Classification hierarchy: MCS > Cloud Types > Unclassified
    # If a cell has an MCS at any time during the window, MCS takes priority
    cloud_types_aggregated = np.where(combined_mcs_swath > 0, 0, cloud_types_aggregated)
    
    # Filter precipitation by cloud type: set to 0 where MCS swath exists
    # This makes the 4 cloud types mutually exclusive with MCS precipitation
    mcs_mask = combined_mcs_swath > 0
    dc_pr = np.where(mcs_mask, 0, dc_pr)
    st_pr = np.where(mcs_mask, 0, st_pr)
    nd_pr = np.where(mcs_mask, 0, nd_pr)
    dz_pr = np.where(mcs_mask, 0, dz_pr)

    # import matplotlib.pyplot as plt
    # cloud_types_da = xr.DataArray(cloud_types_aggregated, dims=['cell'], coords={'lat': ('cell', _ds['lat'].values), 'lon': ('cell', _ds['lon'].values), 'cell': _ds['cell'].values})
    # egh.healpix_show(cloud_types_da.where(cloud_types_da > 0), cmap='tab10', vmin=0, vmax=3)
    # import pdb; pdb.set_trace()
    if verbose:
        n_classified = (cloud_types_aggregated > 0).sum()
        n_mcs_overlap = n_classified_before - n_classified
        print(f"  Cloud types after MCS filter: {n_classified:,} cells")
        print(f"  Removed {n_mcs_overlap:,} cells due to MCS overlap")
        print(f"  ✅ Completed processing for this time chunk")

    return {
        'mcs_mask': combined_mcs_swath,
        'cloud_types': cloud_types_aggregated,
        'dc_pr': dc_pr,
        'st_pr': st_pr,
        'nd_pr': nd_pr,
        'dz_pr': dz_pr,
    }

def process_timechunk_wrapper_zarr(start_idx, end_idx, zarr_path, config, verbose=False,
                                    store_offsets=None):
    """
    Wrapper function for processing a time chunk by reading from zarr file.
    
    This approach avoids serialization issues by:
    1. Having each worker read the zarr file directly
    2. Using integer indices instead of passing time arrays (reduces graph size)
    
    Args:
        start_idx: Start index in the time dimension
        end_idx: End index in the time dimension (exclusive)
        zarr_path: Path (str) or list of paths to the zarr store(s).
            When a list is provided the stores are opened and concatenated
            along the time dimension before slicing.
        config: Configuration dictionary with catalog and variable information
        verbose: Whether to print verbose output
        store_offsets: Optional list of (global_start, global_end_exclusive, path) tuples.
            When provided, each worker opens only the store(s) that contain the
            requested time range instead of opening and concatenating all stores.
        
    Returns:
        tuple: (time_str, results_dict) or (time_str, None) if error
    """
    try:
        # Open only the zarr store(s) needed for this time slice.
        # store_offsets lets each worker target a single annual store rather than
        # opening and concatenating all stores (avoids per-worker memory overhead).
        open_stores = []
        if store_offsets is not None:
            relevant = [(gs, ge, p) for gs, ge, p in store_offsets
                        if gs < end_idx and ge > start_idx]
            if len(relevant) == 1:
                gs, _ge, path = relevant[0]
                _s = xr.open_dataset(path, engine='zarr', chunks=None)
                open_stores.append(_s)
                _ds = _s.isel(time=slice(start_idx - gs, end_idx - gs))
            else:
                parts = []
                for gs, ge, path in relevant:
                    _s = xr.open_dataset(path, engine='zarr', chunks=None)
                    open_stores.append(_s)
                    local_start = max(0, start_idx - gs)
                    local_end = min(ge - gs, end_idx - gs)
                    parts.append(_s.isel(time=slice(local_start, local_end)))
                _ds = xr.concat(parts, dim='time')
        elif isinstance(zarr_path, list):
            for p in zarr_path:
                _s = xr.open_dataset(p, engine='zarr', chunks=None)
                open_stores.append(_s)
            _ds = xr.concat(open_stores, dim='time').isel(time=slice(start_idx, end_idx))
        else:
            _s = xr.open_dataset(zarr_path, engine='zarr', chunks=None)
            open_stores.append(_s)
            _ds = _s.isel(time=slice(start_idx, end_idx))
        
        # Get the first time value for output
        out_time_val = _ds.time.values[0]
        
        # Load zarr data into memory and close the open stores to free memory
        _ds = _ds.load()
        _ds = _ds.pipe(egh.attach_coords)
        for _s in open_stores:
            _s.close()
        
        # Add tb and pr variables to the dataset chunk
        _ds = add_tb_pr_to_dataset(_ds, config)
        # Eagerly load tb/pr into numpy to release the catalog's dask task graphs.
        # Without this, each .values call in process_timechunk_swath triggers a
        # separate catalog read and keeps the full task graph alive in worker memory.
        _ds = _ds.load()
        gc.collect()
        
        # Process this time chunk (all times in the chunk)
        timestep_results = process_timechunk_swath(_ds, verbose=verbose)
        
        # Explicit cleanup
        del _ds
        gc.collect()
        
        return str(out_time_val), timestep_results
        
    except Exception as e:
        print(f"Error processing time chunk {start_idx}-{end_idx}: {e}")
        traceback.print_exc()
        return None, None

#--------------------------------------------------------------------------------------------------
def stream_process_to_zarr(time_coords, mask_variables, output_path,
                          time_groups, output_time_coords,
                          client=None, logger=None, parallel=True, 
                          input_zarr_path=None, batch_size=100,
                          config=None, store_offsets=None):
    """
    Stream process time chunks and write results to zarr with optional parallel processing.
    
    Uses a batched submission approach to avoid overwhelming the Dask scheduler with
    too many tasks at once. Instead of submitting all chunks at once, submits them
    in batches (super-chunks), waits for each batch to complete, then moves to the next.
    
    Parameters:
    -----------
    time_coords : array-like
        Array of input time coordinates to process
    mask_variables : list
        List of mask variable names to create
    output_path : str
        Path to output zarr file
    time_groups : dict
        Dictionary mapping aligned output times to lists of input time indices
    output_time_coords : array-like
        Array of aligned output time coordinates
    client : dask.distributed.Client, optional
        Dask client for parallel processing
    logger : logging.Logger, optional
        Logger instance
    parallel : bool
        Whether to use parallel processing
    input_zarr_path : str
        Path to input zarr file for workers to read from
    batch_size : int
        Number of output chunks to submit per batch (default: 100)
    config : dict
        Configuration dictionary with catalog and variable information
        
    Returns:
    --------
    int : Number of successfully processed chunks
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Time-aligned processing using time groups
    total_chunks = len(output_time_coords)
    output_times_list = list(output_time_coords)
    chunk_indices = list(range(total_chunks))
    
    logger.info(f"Processing {len(time_coords)} input time steps into {total_chunks} aligned time groups")
    logger.info(f"Each chunk aggregates multiple hourly time steps into 1 swath mask at standard hours")
    
    # Process and write time steps in chunks
    total_processed = 0
    
    if parallel and client is not None:
        # PARALLEL MODE: Submit chunks in batches to avoid overwhelming scheduler
        total_batches = (len(chunk_indices) + batch_size - 1) // batch_size
        logger.info(f"Using batched submission: {total_batches} batches of up to {batch_size} chunks each")
        
        # Prepare chunk metadata for chunks we're actually processing
        chunk_metadata = []
        for chunk_idx in chunk_indices:
            # Time-aligned processing: get indices from time groups
            aligned_time = output_times_list[chunk_idx]
            # Convert numpy datetime64 to pandas Timestamp for dictionary lookup
            aligned_time_pd = pd.Timestamp(aligned_time)
            time_indices = time_groups[aligned_time_pd]
            start_idx = int(time_indices[0])
            end_idx = int(time_indices[-1]) + 1
            output_time = aligned_time
            
            chunk_metadata.append({
                'chunk_idx': chunk_idx,
                'start_idx': start_idx,
                'end_idx': end_idx,
                'output_time': output_time
            })
        
        # Process chunks in batches
        for batch_idx in range(total_batches):
            batch_start = batch_idx * batch_size
            batch_end = min((batch_idx + 1) * batch_size, len(chunk_metadata))
            batch_chunks = chunk_metadata[batch_start:batch_end]
            
            logger.info(f"")
            logger.info(f"{'='*80}")
            logger.info(f"BATCH {batch_idx + 1}/{total_batches}: Processing {len(batch_chunks)} chunks")
            logger.info(f"{'='*80}")
            
            # Submit all chunks in this batch to Dask workers (using indices only)
            futures = {}
            for meta in batch_chunks:
                future = client.submit(
                    process_timechunk_wrapper_zarr,
                    meta['start_idx'],
                    meta['end_idx'],
                    input_zarr_path,
                    config,
                    verbose=False,
                    store_offsets=store_offsets
                )
                futures[future] = meta
            
            logger.info(f"Submitted {len(futures)} chunks to workers for this batch...")
            
            # Process results as they complete within this batch
            for future in as_completed(futures):
                meta = futures[future]
                chunk_idx = meta['chunk_idx']
                start_idx = meta['start_idx']
                end_idx = meta['end_idx']
                
                logger.info(f"Processing chunk {chunk_idx + 1}/{total_chunks}: time steps {start_idx}-{end_idx-1}")
                
                chunk_results = {}
                try:
                    time_str, result = future.result()
                    if result is not None and time_str is not None:
                        chunk_results[time_str] = result
                    else:
                        logger.warning(f"Skipping chunk {chunk_idx + 1} (time steps {start_idx}-{end_idx-1}) due to processing error")
                except Exception as e:
                    logger.error(f"Error in Dask task for chunk {chunk_idx + 1}: {e}")
                    traceback.print_exc()
                    continue
                
                # Write this chunk to zarr immediately
                if len(chunk_results) > 0:
                    try:
                        # Get the output time for this chunk
                        # Use aligned time from metadata (already computed)
                        output_time = np.array([meta['output_time']], dtype='datetime64[ns]')
                        
                        logger.info(f"Writing chunk {chunk_idx + 1} to zarr...")
                        append_chunk_to_zarr(
                            chunk_results=chunk_results,
                            chunk_times=output_time,
                            chunk_idx=chunk_idx,
                            mask_variables=mask_variables,
                            output_path=output_path,
                            logger=logger
                        )
                        
                        # Update progress
                        processed_this_chunk = len(chunk_results)
                        total_processed += processed_this_chunk
                        logger.info(f"Chunk {chunk_idx + 1} complete: {processed_this_chunk} swath mask(s) written")
                        
                        # Free memory
                        del chunk_results
                        gc.collect()
                        
                    except Exception as e:
                        logger.error(f"Error writing chunk {chunk_idx + 1} to zarr: {e}")
                        traceback.print_exc()
                        continue
                else:
                    logger.warning(f"No valid results for chunk {chunk_idx + 1}, skipping write")
            
            logger.info(f"Batch {batch_idx + 1}/{total_batches} complete: {total_processed}/{len(chunk_indices)} total chunks processed")
    
    else:
        # SERIAL MODE: Process chunks one at a time
        logger.info("Processing chunks in serial mode...")
        
        for chunk_idx in chunk_indices:
            # Time-aligned processing: get indices from time groups
            aligned_time = output_times_list[chunk_idx]
            # Convert numpy datetime64 to pandas Timestamp for dictionary lookup
            aligned_time_pd = pd.Timestamp(aligned_time)
            time_indices = time_groups[aligned_time_pd]
            start_idx = int(time_indices[0])
            end_idx = int(time_indices[-1]) + 1
            output_time_val = aligned_time
            
            logger.info(f"Processing chunk {chunk_idx + 1}/{total_chunks}: time steps {start_idx}-{end_idx-1}")
            
            chunk_results = {}
            time_str, result = process_timechunk_wrapper_zarr(
                start_idx, end_idx, input_zarr_path, config,
                verbose=False, store_offsets=store_offsets
            )
            if result is not None and time_str is not None:
                chunk_results[time_str] = result
            else:
                logger.warning(f"Skipping chunk {chunk_idx + 1} (time steps {start_idx}-{end_idx-1}) due to processing error")
            
            # Write this chunk to zarr immediately
            if len(chunk_results) > 0:
                try:
                    # Get the output time for this chunk (aligned time)
                    output_time = np.array([output_time_val], dtype='datetime64[ns]')
                    
                    logger.info(f"Writing chunk {chunk_idx + 1} to zarr...")
                    append_chunk_to_zarr(
                        chunk_results=chunk_results,
                        chunk_times=output_time,
                        chunk_idx=chunk_idx,
                        mask_variables=mask_variables,
                        output_path=output_path,
                        logger=logger
                    )
                    
                    # Update progress
                    processed_this_chunk = len(chunk_results)
                    total_processed += processed_this_chunk
                    logger.info(f"Chunk {chunk_idx + 1} complete: {processed_this_chunk} swath mask(s) written")
                    
                    # Free memory
                    del chunk_results
                    gc.collect()
                    
                except Exception as e:
                    logger.error(f"Error writing chunk {chunk_idx + 1} to zarr: {e}")
                    traceback.print_exc()
                    continue
            else:
                logger.warning(f"No valid results for chunk {chunk_idx + 1}, skipping write")

    logger.info(f"Stream processing complete: {total_processed}/{len(chunk_indices)} chunks written successfully")
    return total_processed

#--------------------------------------------------------------------------------------------------
def main():
    """Main function to run the make MCS swath process"""

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process MCS swath masks')
    parser.add_argument("-c", "--config", help="yaml config file for PyFLEXTRKR MCS tracking", required=True)
    parser.add_argument('--parallel', action='store_true', default=True,
                       help='Use parallel processing with Dask (default: True)')
    parser.add_argument('--no-parallel', action='store_false', dest='parallel',
                       help='Disable parallel processing')
    parser.add_argument('--workers', type=int, default=32,
                       help='Number of Dask workers (default: 32)')
    parser.add_argument('--threads-per-worker', type=int, default=1,
                       help='Number of threads per worker (default: 1)')
    parser.add_argument('--aggregation-window', type=int, default=6,
                       help='Number of hourly time steps to aggregate into one swath mask (default: 6)')
    parser.add_argument('--batch-size', type=int, default=100,
                       help='Number of chunks to submit per batch to avoid overwhelming scheduler (default: 100)')
    parser.add_argument('--test-steps', type=int, default=None,
                       help='Number of time steps to process for testing (default: all)')
    
    args = parser.parse_args()

    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)

    start_time = time.time()
    logger.info("Starting making MCS mask swath...")

    # Get aggregation window from input argument 
    # This determines the time window size to aggregate MCS masks (e.g., 6 = aggregate 6 hourly steps into 1 swath)
    aggregation_window = args.aggregation_window
    
    # Zarr storage chunk size (number of output time steps per zarr chunk)
    # 28 = 1 week of 6-hourly data (7 days * 4 swaths/day)
    zarr_chunk_size_time = 28
    
    # Define output variables
    mask_variables = ['mcs_mask', 'cloud_types', 'dc_pr', 'st_pr', 'nd_pr', 'dz_pr']
    
    # Define attributes for each variable
    var_attrs = {
        'mcs_mask': {
            'long_name': 'MCS swath mask',
            'description': 'MCS track swath mask aggregated over time window. Each pixel labeled with MCS track number.',
            'units': '1',
            'valid_range': [0, 1000000],
            'comment': 'Track numbers are assigned by PyFLEXTRKR. Overlapping tracks resolved by highest coverage count.'
        },
        'cloud_types': {
            'long_name': 'Cloud type classification',
            'description': 'Priority-based cloud type within non-MCS cold cloud shield over time window',
            'units': '1',
            'valid_range': [0, 4],
            'flag_values': [0, 1, 2, 3, 4],
            'flag_meanings': 'no_cloud deep_convective stratiform non_deep_convective drizzle',
            'comment': 'Classification based on brightness temperature and precipitation rate. Priority ranking: 1>2>3>4. Mutually exclusive with MCS.'
        },
        'dc_pr': {
            'long_name': 'Deep convective precipitation',
            'description': 'Frequency-weighted mean precipitation for deep convective clouds',
            'units': 'mm h-1',
            'comment': 'Mean precipitation when cell is classified as deep convective (type 1), weighted by frequency. Mutually exclusive with MCS swath.'
        },
        'st_pr': {
            'long_name': 'Stratiform precipitation',
            'description': 'Frequency-weighted mean precipitation for stratiform clouds',
            'units': 'mm h-1',
            'comment': 'Mean precipitation when cell is classified as stratiform (type 2), weighted by frequency. Mutually exclusive with MCS swath.'
        },
        'nd_pr': {
            'long_name': 'Non-deep convective precipitation',
            'description': 'Frequency-weighted mean precipitation for non-deep convective clouds',
            'units': 'mm h-1',
            'comment': 'Mean precipitation when cell is classified as non-deep convective (type 3), weighted by frequency. Mutually exclusive with MCS swath.'
        },
        'dz_pr': {
            'long_name': 'Drizzle precipitation',
            'description': 'Frequency-weighted mean precipitation for drizzle',
            'units': 'mm h-1',
            'comment': 'Mean precipitation when cell is classified as drizzle (type 4), weighted by frequency. Mutually exclusive with MCS swath.'
        }
    }
    
    # Parallel processing configuration
    parallel = args.parallel
    n_workers = args.workers
    threads_per_worker = args.threads_per_worker
    batch_size = args.batch_size
    config_file = args.config

    # Load MCS tracking configuration for the specified source
    zoom = 8
    config = load_config(config_file)
    root_path = config.get("root_path")
    pixel_path_name = config.get("pixel_path_name")
    pixel_path = f"{root_path}{pixel_path_name}/"
    in_basename = config.get("zarr_output_presets", {}).get("healpix").get("out_filebase")
    # # Prefer dated zarr stores (e.g., *_hp8_v1_20190101.0000_20200101.0100.zarr).
    # # Pattern matches *_v1_{date}_*.zarr but not *_v1.zarr.
    # # Falls back to the canonical *_hp8_v1.zarr when no dated stores exist.
    # zarr_pattern = f"{pixel_path}{in_basename}hp{zoom}_v1_[0-9]*.zarr"
    # dated_zarr_stores = sorted(glob.glob(zarr_pattern))
    # if len(dated_zarr_stores) >= 1:
    #     logger.info(f"Found {len(dated_zarr_stores)} dated zarr store(s) to use:")
    #     for p in dated_zarr_stores:
    #         logger.info(f"  {os.path.basename(p)}")
    #     in_zarr = dated_zarr_stores if len(dated_zarr_stores) > 1 else dated_zarr_stores[0]
    # else:
    # Fall back to the canonical v1 store
    in_zarr = f"{pixel_path}{in_basename}hp{zoom}_v1.zarr"
    logger.info(f"No dated zarr stores found, using: {os.path.basename(in_zarr)}")

    # Catalog information
    catalog_file = config.get('catalog_file')
    catalog_location = config.get('catalog_location')
    catalog_source = config.get('catalog_source')
    catalog_params = config.get('catalog_params', {})
    # Update the zoom in the catalog_params to match the zoom in the mask file
    catalog_params.update({'zoom': zoom})
    # Precipitation, OLR variable names and conversion factor
    varname_precip_liq = 'pr'
    varname_precip_ice = 'prs'
    varname_olr = 'rlut'
    pr_convert_factor = config.get('pcp_convert_factor')

    # Get source name from root path (e.g., /pscratch/sd/w/wcmca1/hackathon/mcs/scream/)
    source_name = os.path.basename(os.path.normpath(root_path))
    # Strip trailing year-range suffix (e.g., IMERGv7_2019_2021 -> IMERGv7)
    source_name = re.sub(r'_20\d{2}_20\d{2}$', '', source_name)

    # Output paths
    out_dir = "/pscratch/sd/w/wcmca1/hackathon/mcs_masks/"
    out_basename = f"{source_name}_mcs_masks_hp{zoom}.zarr"
    out_zarr = f"{out_dir}{out_basename}"
    os.makedirs(out_dir, exist_ok=True)

    # Prepare config dictionary for processing
    processing_config = {
        'catalog_file': catalog_file,
        'catalog_location': catalog_location,
        'catalog_source': catalog_source,
        'catalog_params': catalog_params,
        'varname_olr': varname_olr,
        'varname_precip_liq': varname_precip_liq,
        'varname_precip_ice': varname_precip_ice,
        'pcp_convert_factor': pr_convert_factor,
    }

    print("="*80)
    print("MAKE MCS SWATH MASK PROCESSING")
    print("="*80)
    print(f"Source: {source_name}")
    if isinstance(in_zarr, list):
        print(f"Input: {len(in_zarr)} zarr stores (concatenated by time):")
        for p in in_zarr:
            print(f"  {p}")
    else:
        print(f"Input: {in_zarr}")
    print(f"Output: {out_zarr}")
    print(f"Parallel processing: {parallel}")
    if parallel:
        print(f"Workers: {n_workers}, Threads per worker: {threads_per_worker}")
        print(f"Batch size: {batch_size} chunks per batch")

    # Setup Dask client
    client = setup_dask_client(
        parallel=parallel, 
        n_workers=n_workers, 
        threads_per_worker=threads_per_worker, 
        logger=logger,
    )

    try:
        # Load the full dataset
        print(f"\nLoading full dataset...")
        store_offsets = None  # per-store offset table so workers open only what they need
        try:
            # Read MCS mask dataset from zarr (single store or multiple dated stores).
            # Build store_offsets alongside so workers can target individual stores.
            if isinstance(in_zarr, list):
                ds_list = []
                store_offsets = []
                offset = 0
                for p in in_zarr:
                    _s = xr.open_zarr(p, consolidated=True, mask_and_scale=True)
                    n = len(_s.time)
                    store_offsets.append((offset, offset + n, p))
                    offset += n
                    ds_list.append(_s)
                ds = xr.concat(ds_list, dim='time')
            else:
                ds = xr.open_zarr(in_zarr, consolidated=True, mask_and_scale=True)
            
            # ds = ds.pipe(egh.attach_coords)  # Commented out for testing
            print(f"  ✅ Dataset loaded successfully")
            print(f"  Time steps: {len(ds.time)}")
            print(f"  Data variables: {list(ds.data_vars)}")
            print(f"  Spatial dimensions: {dict(ds.dims)}")
        
        except Exception as e:
            print(f"  ❌ Error loading dataset: {e}")
            return
        # Limit time steps for testing if requested
        if args.test_steps is not None:
            ds = ds.isel(time=slice(0, args.test_steps))
            print(f"  📋 Limited to {args.test_steps} time steps for testing")
            
        # Initialize streaming processing configuration
        print(f"\nInitializing streaming zarr processing...")
        time_coords = ds["time"].values

        # Add processing metadata
        attrs = ds.attrs.copy()
        attrs.update({
            'processing_info': 'MCS swath mask creation',
            'processing_date': str(np.datetime64('today')),
            'processing_script': os.path.basename(__file__),
        })

        # Stream processing and writing to zarr
        print(f"\nStreaming processing and writing {len(time_coords)} time steps to zarr...")
        
        print(f"Using aggregation_window={aggregation_window} hours for swath computation")

        # Create output time coordinates aligned to standard hours (00, 06, 12, 18)
        # Convert to pandas for easier time manipulation
        time_df = pd.DataFrame({'time': pd.to_datetime(time_coords)})
        
        # Round times to nearest aggregation_window hour boundary
        # For 6-hour aggregation: rounds to 00, 06, 12, 18
        time_df['time_aligned'] = time_df['time'].dt.floor(f'{aggregation_window}h')
        
        # Get unique aligned times (these are the output time coordinates)
        output_time_coords = time_df['time_aligned'].unique()
        output_time_coords = np.array(output_time_coords, dtype='datetime64[ns]')
        
        # Create mapping from aligned times to input time indices
        # This tells us which input times belong to each output time
        time_groups = time_df.groupby('time_aligned').groups
        
        logger.info(f"Input: {len(time_coords)} hourly time steps")
        logger.info(f"Output: {len(output_time_coords)} {aggregation_window}-hourly swath times (aligned to standard hours)")
        logger.info(f"First input time: {time_coords[0]}, First output time: {output_time_coords[0]}")
        logger.info(f"Last input time: {time_coords[-1]}, Last output time: {output_time_coords[-1]}")
        logger.info(f"Zarr time chunking: {zarr_chunk_size_time} time steps per chunk (optimized for weekly access)")

        try:
            # Initialize the zarr store structure (once only)
            logger.info("Initializing zarr store...")
            initialize_zarr_store(
                output_path=out_zarr,
                time_coords=output_time_coords,  # Use aggregated time coordinates
                mask_variables=mask_variables,
                template_coords=ds.coords,
                attrs=attrs,
                chunk_size_time=zarr_chunk_size_time,  # Zarr storage chunks (28 = 1 week)
                var_attrs=var_attrs  # Variable attributes
            )

            # Stream process with chunked zarr writing
            total_processed = stream_process_to_zarr(
                time_coords=time_coords,  # Pass input time coords for processing
                mask_variables=mask_variables,
                output_path=out_zarr,
                time_groups=time_groups,  # Mapping of aligned times to input indices
                output_time_coords=output_time_coords,  # Aligned output time coordinates
                client=client,
                logger=logger,
                parallel=parallel,
                input_zarr_path=in_zarr,  # Pass the input zarr path for workers
                batch_size=batch_size,  # Number of chunks per batch
                config=processing_config,  # Configuration dictionary
                store_offsets=store_offsets,  # Per-store offset table for targeted reads
            )
            
            logger.info(f"✅ Processing complete: {total_processed} chunks written to {out_zarr}")

        except Exception as e:
            logger.error(f"Error writing chunked zarr: {e}")
            traceback.print_exc()
            print(f"  ❌ Error writing zarr: {e}")
            return

        # Calculate total processing time
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        # Store success info to print after Dask cleanup
        success_info = {
            'elapsed_time': elapsed_time,
            'out_zarr': out_zarr,
            'total_processed': total_processed
        }


    finally:
        # Always cleanup client
        if client and parallel:
            # Suppress Dask shutdown messages by temporarily raising log level
            logging.getLogger('distributed').setLevel(logging.CRITICAL)
            logging.getLogger('distributed.worker').setLevel(logging.CRITICAL)
            logging.getLogger('distributed.nanny').setLevel(logging.CRITICAL)
            
            client.close()
        
        # Print success message after Dask cleanup (so it's always visible at the end)
        if 'success_info' in locals():
            elapsed = success_info['elapsed_time']
            print(f"\n{'='*80}")
            print(f"✅ PROCESSING COMPLETE!")
            print(f"{'='*80}")
            print(f"Output: {success_info['out_zarr']}")
            print(f"Chunks processed: {success_info['total_processed']}")
            print(f"Total time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
            print(f"{'='*80}\n")

if __name__ == "__main__":
    main()