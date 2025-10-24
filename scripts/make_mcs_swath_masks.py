import xarray as xr
import numpy as np
import pandas as pd
import cftime
import yaml
import calendar
import os, glob
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

#-------------------------------------------------------------------
def convert_cftime_to_standard_calendar(cftime_times):
    """
    Convert cftime datetime objects (non-standard calendars) to numpy datetime64 (standard calendar).
    
    This function handles conversion from calendars like DatetimeNoLeap, Datetime360Day, etc.
    to standard proleptic_gregorian calendar (numpy.datetime64).
    
    Args:
        cftime_times: array-like of cftime datetime objects
            Timestamps with non-standard calendar (e.g., DatetimeNoLeap, Datetime360Day)
            
    Returns:
        numpy.ndarray: Array of numpy.datetime64 objects with standard calendar
    """
    
    # Check if input is a single timestamp
    is_single_object = not hasattr(cftime_times, '__iter__')
    
    # Convert to list for uniform processing
    times_list = [cftime_times] if is_single_object else cftime_times
    
    # Convert each cftime timestamp to numpy datetime64
    converted_times = []
    for t in times_list:
        # Check if it's already a standard datetime type
        if isinstance(t, (np.datetime64, pd.Timestamp)):
            converted_times.append(np.datetime64(t, 'ns'))
        # If it's a cftime object, extract components and create datetime64
        elif hasattr(t, 'year'):
            # Create a pandas Timestamp from the cftime components
            # Note: This may shift dates for non-standard calendars that have different
            # day counts (e.g., Feb 30 in 360_day calendar doesn't exist in standard)
            try:
                pd_time = pd.Timestamp(
                    year=t.year, month=t.month, day=t.day,
                    hour=t.hour, minute=t.minute, second=t.second
                )
                converted_times.append(pd_time.to_datetime64())
            except ValueError as e:
                # Handle invalid dates (e.g., Feb 30)
                # For simplicity, we'll skip invalid dates or adjust them
                print(f"Warning: Could not convert {t} to standard calendar: {e}")
                # Try to adjust the day to the last valid day of the month
                last_day = calendar.monthrange(t.year, t.month)[1]
                adjusted_day = min(t.day, last_day)
                pd_time = pd.Timestamp(
                    year=t.year, month=t.month, day=adjusted_day,
                    hour=t.hour, minute=t.minute, second=t.second
                )
                converted_times.append(pd_time.to_datetime64())
        else:
            raise TypeError(f"Unsupported time type: {type(t)}")
    
    # Return a single object or array based on input type
    if is_single_object:
        return converted_times[0]
    else:
        return np.array(converted_times, dtype='datetime64[ns]')


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
def classify_cloud_types(tb, pr, tb_thresh, nonmcs_ccs_mask, pr_threshold=0.5):
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
    nonmcs_ccs_mask : numpy.ndarray
        Non-MCS CCS mask (only classify where > 0)
    pr_threshold : float, optional
        Precipitation threshold for classification (default: 0.5 mm/h)
    
    Returns:
    --------
    cloud_type : numpy.ndarray
        Cloud type classification:
        0 = No cloud (or outside nonmcs_ccs_mask)
        1 = Deep convective (tb < tb_thresh & pr >= pr_threshold)
        2 = Stratiform (tb < tb_thresh & pr < pr_threshold)
        3 = Non-deep convective (tb >= tb_thresh & pr >= pr_threshold)
        4 = Drizzle (tb >= tb_thresh & pr < pr_threshold)
    """
    # Initialize cloud type array with zeros
    cloud_type = np.zeros_like(tb, dtype=np.int8)
    
    # Only classify where nonmcs_ccs_mask > 0
    valid_mask = nonmcs_ccs_mask > 0
    
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
def find_most_frequent_cloud_type(cloud_types_time_series):
    """
    Find the most frequent cloud type over time for each cell (VECTORIZED).
    
    For ties, prioritize lower type numbers (e.g., 1 > 2 > 3 > 4).
    Automatically determines the number of cloud types from the data.
    
    This vectorized version uses pure numpy operations to process all cells 
    simultaneously, providing dramatic speedup (~135x faster) compared to 
    loop-based approaches. Essential for processing millions of HEALPix cells.
    
    Parameters:
    -----------
    cloud_types_time_series : numpy.ndarray
        Array of shape (n_times, n_cells) with cloud type values starting from 0
        Type 0 is assumed to be "no cloud" and is excluded from results
        
    Returns:
    --------
    most_frequent : numpy.ndarray
        Array of shape (n_cells,) with the most frequent cloud type
    """
    n_times, n_cells = cloud_types_time_series.shape
    
    # Determine the number of unique cloud types from the data
    # max_type = int(cloud_types_time_series.max())
    max_type = 4
    n_types = max_type + 1  # e.g., if max is 4, we have types 0-4 (5 types)
    
    # Count occurrences of each type for all cells at once
    # Initialize count array: shape (n_types, n_cells)
    counts = np.zeros((n_types, n_cells), dtype=np.int16)
    
    # Count occurrences for each type using vectorized operations
    for type_val in range(n_types):
        counts[type_val, :] = (cloud_types_time_series == type_val).sum(axis=0)
    
    # For each cell, find the type with maximum count (excluding type 0)
    # counts[1:, :] excludes type 0
    cloud_counts = counts[1:, :]
    
    # Find the maximum count for each cell (across all cloud types except 0)
    max_counts = cloud_counts.max(axis=0)  # shape: (n_cells,)
    
    # Initialize result array
    most_frequent = np.zeros(n_cells, dtype=np.int8)
    
    # For cells with clouds, find the type with max count
    # In case of ties, select the smallest type number (highest priority)
    cells_with_clouds = max_counts > 0
    
    if cells_with_clouds.any():
        # For each cell with clouds, find which type has the max count
        # Check in priority order (1, 2, 3, ...) to handle ties
        for type_val in range(1, n_types):
            # Cells where this type has the max count
            is_max = (cloud_counts[type_val - 1, :] == max_counts) & cells_with_clouds
            # Set these cells to this type (only if not already set)
            most_frequent[is_max & (most_frequent == 0)] = type_val
    
    return most_frequent

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
    
    # 3. Select ds_p to only the common times (BEFORE expensive operations)
    ds_p = ds_p.sel(time=common_times)
    
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
        ccs_mask_t = _ds['ccs_mask'].isel(time=t).values
        
        # Replace NaN with 0 (when mask_and_scale=True, fill_value becomes NaN)
        mcs_mask_t = np.nan_to_num(mcs_mask_t, nan=0.0).astype(int)
        ccs_mask_t = np.nan_to_num(ccs_mask_t, nan=0.0).astype(int)
        
        # Compute non-MCS CCS mask for this time step
        nmcs_ccs_mask_t = np.where(mcs_mask_t == 0, ccs_mask_t, 0)
        
        # Classify cloud types
        cloud_type_t = classify_cloud_types(tb_t, pr_t, tb_thresh, nmcs_ccs_mask_t, pr_threshold=0.5)
        cloud_types_timeseries.append(cloud_type_t)
    
    # Stack into array (time, cell)
    cloud_types_timeseries = np.stack(cloud_types_timeseries, axis=0)
    
    # Find the most frequent cloud type over the time window
    cloud_types_aggregated = find_most_frequent_cloud_type(cloud_types_timeseries)

    # import matplotlib.pyplot as plt
    # cloud_types_da = xr.DataArray(cloud_types_aggregated, dims=['cell'], coords={'lat': ('cell', _ds['lat'].values), 'lon': ('cell', _ds['lon'].values), 'cell': _ds['cell'].values})
    # egh.healpix_show(cloud_types_da.where(cloud_types_da > 0), cmap='tab10', vmin=0, vmax=3)
    # import pdb; pdb.set_trace()
    if verbose:
        n_classified = (cloud_types_aggregated > 0).sum()
        print(f"  ✅ Completed processing for this time chunk")
        print(f"     Classified {n_classified:,} cells with cloud types")

    return {
        'mcs_mask': combined_mcs_swath,
        # 'ccs_mask': ccs_mask_sum,
        'cloud_types': cloud_types_aggregated,
    }

def process_timechunk_wrapper_zarr(start_idx, end_idx, zarr_path, config, verbose=False):
    """
    Wrapper function for processing a time chunk by reading from zarr file.
    
    This approach avoids serialization issues by:
    1. Having each worker read the zarr file directly
    2. Using integer indices instead of passing time arrays (reduces graph size)
    
    Args:
        start_idx: Start index in the time dimension
        end_idx: End index in the time dimension (exclusive)
        zarr_path: Path to the zarr file containing the data
        config: Configuration dictionary with catalog and variable information
        verbose: Whether to print verbose output
        
    Returns:
        tuple: (time_str, results_dict) or (time_str, None) if error
    """
    try:
        # Each worker opens the zarr file independently
        ds = xr.open_dataset(zarr_path, engine='zarr', chunks=None)
        
        # Select the specific times for this chunk using integer indexing
        _ds = ds.isel(time=slice(start_idx, end_idx))
        
        # Get the first time value for output
        out_time_val = _ds.time.values[0]
        
        # Load the data into memory
        _ds = _ds.load()
        _ds = _ds.pipe(egh.attach_coords)
        
        # Close the full dataset to free memory
        ds.close()
        
        # Add tb and pr variables to the dataset chunk
        _ds = add_tb_pr_to_dataset(_ds, config)
        
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
                          config=None):
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
                    verbose=False
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
                verbose=False
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
    mask_variables = ['mcs_mask', 'cloud_types']
    
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
            'description': 'Most frequent cloud type within non-MCS cold cloud shield over time window',
            'units': '1',
            'valid_range': [0, 4],
            'flag_values': [0, 1, 2, 3, 4],
            'flag_meanings': 'no_cloud deep_convective stratiform non_deep_convective drizzle',
            'comment': 'Classification based on brightness temperature and precipitation rate. Priority for ties: 1>2>3>4.'
        }
    }
    
    # Parallel processing configuration
    parallel = args.parallel
    n_workers = args.workers
    threads_per_worker = args.threads_per_worker
    batch_size = args.batch_size
    config_file = args.config

    # Load configuration for the specified source
    zoom = 8
    config = load_config(config_file)
    root_path = config.get("root_path")
    pixel_path_name = config.get("pixel_path_name")
    pixel_path = f"{root_path}{pixel_path_name}/"
    in_basename = config.get("zarr_output_presets", {}).get("healpix").get("out_filebase")
    in_zarr = f"{pixel_path}{in_basename}hp{zoom}_v1.zarr"
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
        try:
            # Read MCS mask dataset from zarr
            # ds = xr.open_zarr(in_zarr, consolidated=True, mask_and_scale=False)
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