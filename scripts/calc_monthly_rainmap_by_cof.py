import numpy as np
import sys, os
from pathlib import Path
# Add src directory to path for zarr_tools import
sys.path.append(str(Path(__file__).parent.parent / 'src'))
from zarr_tools import setup_dask_client
import yaml
import xarray as xr
import pandas as pd
import time
import psutil
import argparse
import cftime
import intake
import gc
import logging
import easygems.healpix as egh

def parse_cmd_args():
    # Define and retrieve the command-line arguments...
    parser = argparse.ArgumentParser(
        description="Calculate monthly MCS precipitation statistics."
    )
    parser.add_argument("-c", "--config", help="yaml config file for datasets", required=True)
    parser.add_argument("--source", help="catalog source name from config file", required=True)
    # parser.add_argument("--zoom", help="HEALPix zoom level", type=int, default=None)
    parser.add_argument("--workers", help="number of Dask workers (default: 14)", type=int, default=14)
    parser.add_argument("--threads-per-worker", help="threads per worker (default: 4)", type=int, default=4)
    parser.add_argument("--chunk_days", help="number of days to process in each chunk (default: 6)", type=int, default=6)
    parser.add_argument("--pcp_thresh", help="precipitation threshold in mm/h (default: 0.1)", type=float, default=0.1)
    args = parser.parse_args()

    # Put arguments in a dictionary
    args_dict = {
        'config_file': args.config,
        'source': args.source,
        # 'zoom': args.zoom,
        'n_workers': args.workers,
        'threads_per_worker': args.threads_per_worker,
        'chunk_days': args.chunk_days,
        'pcp_thresh': args.pcp_thresh,
    }

    return args_dict

def get_memory_usage():
    """Get current memory usage in a human-readable format"""
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    # Convert to GB
    memory_gb = memory_info.rss / (1024 ** 3)
    return memory_gb

def format_time(seconds):
    """Format time in seconds to hours, minutes, seconds"""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    seconds = int(seconds % 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"

def setup_logging():
    """
    Set the logging message level

    Args:
        None.

    Returns:
        None.
    """
    logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO)


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

def create_union_mask(mask_list):
    """
    Create a union mask from multiple mask variables.
    
    The union operation is applied on the spatial (cell) dimension only, 
    preserving any time or other dimensions in the input masks.
    
    Parameters:
    -----------
    mask_list : list of xarray.DataArray
        List of mask DataArrays to combine. Each mask can contain track IDs (integers > 0)
        or binary values (0 or 1). Masks will be converted to binary before combining.
        All masks should have the same dimensions (e.g., time and cell).
    
    Returns:
    --------
    xarray.DataArray
        Binary mask (0 or 1 as float) representing the union of all input masks.
        A pixel is 1 if it appears in any of the input masks at that time, 0 otherwise.
        The output preserves all dimensions from the input (e.g., time, cell).
    
    Example:
    --------
    >>> # Create union of 3-way overlap masks (with time dimension)
    >>> union_mask = create_union_mask([
    ...     ds.mcs_ar_etc_overlap_mask,  # shape: (time: 3138, cell: 786432)
    ...     ds.ar_mcs_etc_overlap_mask,
    ...     ds.etc_mcs_ar_overlap_mask
    ... ])
    >>> # Result shape: (time: 3138, cell: 786432)
    
    >>> # Create union for single timestep
    >>> union_mask_single = create_union_mask([
    ...     ds.mcs_ar_overlap_mask.sel(time='2020-01-01'),
    ...     ds.ar_etc_overlap_mask.sel(time='2020-01-01')
    ... ])
    
    >>> # Create union of all 2-way overlaps over time range
    >>> all_2way = create_union_mask([
    ...     ds.mcs_ar_overlap_mask.sel(time=slice('2020-01', '2020-02')),
    ...     ds.ar_etc_overlap_mask.sel(time=slice('2020-01', '2020-02')),
    ...     ds.mcs_etc_overlap_mask.sel(time=slice('2020-01', '2020-02'))
    ... ])
    """
    if not mask_list:
        raise ValueError("mask_list cannot be empty")
    
    # Convert each mask to binary (any value > 0 becomes 1)
    # This operation preserves all dimensions (time, cell, etc.)
    binary_masks = [(mask > 0).astype(float) for mask in mask_list]
    
    # Sum all binary masks element-wise across all dimensions
    # This maintains the shape (e.g., time x cell)
    summed_mask = sum(binary_masks)
    
    # Convert back to binary (any sum > 0 becomes 1)
    # Union at each time step: a cell is 1 if it appears in ANY mask at that time
    union_mask = (summed_mask > 0).astype(float)
    
    return union_mask

print("✅ create_union_mask() function defined")

def process_month_chunked(month_ds, chunk_days=5, pcp_thresh=0.1):
    """Process one month of data in time chunks to reduce memory pressure"""
    # Current month's time value for output
    out_time = month_ds.time[0]
    out_time_str = out_time.dt.strftime('%Y-%m').item()
    _year = out_time.dt.year.item()
    _month = out_time.dt.month.item()
    _day = out_time.dt.day.item()
    # Create standard datetime using pandas
    std_time = pd.Timestamp(_year, _month, _day, 0, 0, 0)
    
    # Get total times in month
    ntimes = len(month_ds.time)
    
    # Detect time interval in hours
    if ntimes > 1:
        # Get the first two timestamps
        time_values = month_ds.time.values
        
        # Calculate time interval in hours based on type
        if isinstance(time_values[0], (cftime._cftime.DatetimeNoLeap, cftime._cftime.Datetime360Day)):
            # For cftime objects
            time_interval = (time_values[1] - time_values[0]).total_seconds() / 3600.0
        else:
            # For numpy datetime64 objects
            time_interval = (pd.Timestamp(time_values[1]) - pd.Timestamp(time_values[0])).total_seconds() / 3600.0
            
        print(f"  Detected time interval: {time_interval:.1f} hours")
    else:
        # Default to 1 hours if only one timestamp
        time_interval = 1.0
        print(f"  Using default time interval: {time_interval:.1f} hours")
    
    # Calculate steps per day based on time interval
    steps_per_day = 24.0 / time_interval
    
    # Initialize accumulators for summations
    # Total precipitation
    totprecip_sum = None

    mcs_pcp_sum = None
    
    # Isolated features
    # MCS statistics
    mcs_iso_pcp_sum = None
    mcs_iso_count_sum = None
    mcs_iso_pcpcount_sum = None
    
    # AR statistics
    ar_iso_pcp_sum = None
    ar_iso_count_sum = None
    ar_iso_pcpcount_sum = None
    
    # TC statistics
    tc_pcp_sum = None
    tc_count_sum = None
    tc_pcpcount_sum = None
    
    # ETC statistics
    etc_iso_pcp_sum = None
    etc_iso_count_sum = None
    etc_iso_pcpcount_sum = None

    # Co-occurrence features (2-way)
    mcs_ar_pcp_sum = None
    mcs_ar_count_sum = None
    mcs_ar_pcp_count_sum = None

    mcs_etc_pcp_sum = None
    mcs_etc_count_sum = None
    mcs_etc_pcp_count_sum = None

    ar_etc_pcp_sum = None
    ar_etc_count_sum = None
    ar_etc_pcp_count_sum = None

    # Co-occurrence features (3-way)
    mcs_ar_etc_pcp_sum = None
    mcs_ar_etc_count_sum = None
    mcs_ar_etc_pcp_count_sum = None

    # Cloud types
    # Deep convection (dc)
    dc_pcp_sum = None
    dc_count_sum = None
    
    # Stratiform (st)
    st_pcp_sum = None
    st_count_sum = None
    
    # Non-deep convective (nd)
    nd_pcp_sum = None
    nd_count_sum = None
    
    # Drizzle (dz)
    dz_pcp_sum = None
    dz_count_sum = None

    # Process in chunks of days to limit memory usage
    step = int(chunk_days * steps_per_day)
    for start_idx in range(0, ntimes, step):
        end_idx = min(start_idx + step, ntimes)
        start_day = int(start_idx / steps_per_day) + 1
        end_day = int(end_idx / steps_per_day)
        print(f"  Processing days {start_day}-{end_day} of month {out_time_str}")
        
        # Extract chunk of data and compute immediately to free memory
        chunk_ds = month_ds.isel(time=slice(start_idx, end_idx)).compute()

        # Original masks
        mcs_mask = chunk_ds['mcs_mask']
        ar_mask = chunk_ds['ar_mask']
        etc_mask = chunk_ds['etc_mask']
        tc_mask = chunk_ds['tc_mask']
        cloud_types = chunk_ds['cloud_types']

        # Isolated masks
        mcs_isolated_mask = chunk_ds['mcs_isolated_mask']
        ar_isolated_mask = chunk_ds['ar_isolated_mask']
        etc_isolated_mask = chunk_ds['etc_isolated_mask']

        # Co-occurrence masks
        # Create union of 2-way overlap masks
        mcs_ar_2way_mask = create_union_mask([
            chunk_ds['mcs_ar_overlap_mask'],
            chunk_ds['ar_mcs_overlap_mask'],
        ])
        mcs_etc_2way_mask = create_union_mask([
            chunk_ds['mcs_etc_overlap_mask'],
            chunk_ds['etc_mcs_overlap_mask'],
        ])
        ar_etc_2way_mask = create_union_mask([
            chunk_ds['ar_etc_overlap_mask'],
            chunk_ds['etc_ar_overlap_mask'],
        ])
        # Create union of 3-way overlap masks
        mcs_ar_etc_3way_mask = create_union_mask([
            chunk_ds['mcs_ar_etc_overlap_mask'],
            chunk_ds['ar_mcs_etc_overlap_mask'],
            chunk_ds['etc_mcs_ar_overlap_mask'],
        ])

        precipitation = chunk_ds['pr']
        
        # Compute total precipitation - multiply by time interval to get mm
        chunk_totprecip = (precipitation * time_interval).sum(dim='time')

        # Original features
        chunk_mcs_pcp_sum = (precipitation.where(mcs_mask > 0) * time_interval).sum(dim='time')
        chunk_mcs_count_sum = (mcs_mask > 0).sum(dim='time')
        chunk_mcs_pcp_count_sum = (precipitation.where(mcs_mask > 0) > pcp_thresh).sum(dim='time')

        chunk_ar_pcp_sum = (precipitation.where(ar_mask > 0) * time_interval).sum(dim='time')
        chunk_ar_count_sum = (ar_mask > 0).sum(dim='time')
        chunk_ar_pcp_count_sum = (precipitation.where(ar_mask > 0) > pcp_thresh).sum(dim='time')

        chunk_etc_pcp_sum = (precipitation.where(etc_mask > 0) * time_interval).sum(dim='time')
        chunk_etc_count_sum = (etc_mask > 0).sum(dim='time')
        chunk_etc_pcp_count_sum = (precipitation.where(etc_mask > 0) > pcp_thresh).sum(dim='time')

        # Co-occurrence features
        chunk_mcs_ar_pcp_sum = (precipitation.where(mcs_ar_2way_mask > 0) * time_interval).sum(dim='time')
        chunk_mcs_ar_count_sum = (mcs_ar_2way_mask > 0).sum(dim='time')
        chunk_mcs_ar_pcp_count_sum = (precipitation.where(mcs_ar_2way_mask > 0) > pcp_thresh).sum(dim='time')

        chunk_mcs_etc_pcp_sum = (precipitation.where(mcs_etc_2way_mask > 0) * time_interval).sum(dim='time')
        chunk_mcs_etc_count_sum = (mcs_etc_2way_mask > 0).sum(dim='time')
        chunk_mcs_etc_pcp_count_sum = (precipitation.where(mcs_etc_2way_mask > 0) > pcp_thresh).sum(dim='time')

        chunk_ar_etc_pcp_sum = (precipitation.where(ar_etc_2way_mask > 0) * time_interval).sum(dim='time')
        chunk_ar_etc_count_sum = (ar_etc_2way_mask > 0).sum(dim='time')
        chunk_ar_etc_pcp_count_sum = (precipitation.where(ar_etc_2way_mask > 0) > pcp_thresh).sum(dim='time')

        chunk_mcs_ar_etc_pcp_sum = (precipitation.where(mcs_ar_etc_3way_mask > 0) * time_interval).sum(dim='time')
        chunk_mcs_ar_etc_count_sum = (mcs_ar_etc_3way_mask > 0).sum(dim='time')
        chunk_mcs_ar_etc_pcp_count_sum = (precipitation.where(mcs_ar_etc_3way_mask > 0) > pcp_thresh).sum(dim='time')

        # Isolated features
        # Compute statistics for MCS - multiply by time interval to get mm
        chunk_mcs_iso_pcp_sum = (precipitation.where(mcs_isolated_mask > 0) * time_interval).sum(dim='time')
        chunk_mcs_iso_count_sum = (mcs_isolated_mask > 0).sum(dim='time')
        chunk_mcs_iso_pcp_count_sum = (precipitation.where(mcs_isolated_mask > 0) > pcp_thresh).sum(dim='time')
        
        # Compute statistics for AR - multiply by time interval to get mm
        chunk_ar_iso_pcp_sum = (precipitation.where(ar_isolated_mask > 0) * time_interval).sum(dim='time')
        chunk_ar_iso_count_sum = (ar_isolated_mask > 0).sum(dim='time')
        chunk_ar_iso_pcp_count_sum = (precipitation.where(ar_isolated_mask > 0) > pcp_thresh).sum(dim='time')
        
        # Compute statistics for ETC - multiply by time interval to get mm
        chunk_etc_iso_pcp_sum = (precipitation.where(etc_isolated_mask > 0) * time_interval).sum(dim='time')
        chunk_etc_iso_count_sum = (etc_isolated_mask > 0).sum(dim='time')
        chunk_etc_iso_pcp_count_sum = (precipitation.where(etc_isolated_mask > 0) > pcp_thresh).sum(dim='time')

        # Compute statistics for TC - multiply by time interval to get mm
        chunk_tc_pcp_sum = (precipitation.where(tc_mask > 0) * time_interval).sum(dim='time')
        chunk_tc_count_sum = (tc_mask > 0).sum(dim='time')
        chunk_tc_pcp_count_sum = (precipitation.where(tc_mask > 0) > pcp_thresh).sum(dim='time')
        
        # Cloud types
        # Deep convection (dc): cloud_types == 1
        chunk_dc_pcp_sum = (precipitation.where(cloud_types == 1) * time_interval).sum(dim='time')
        chunk_dc_count_sum = (cloud_types == 1).sum(dim='time')
        
        # Stratiform (st): cloud_types == 2
        chunk_st_pcp_sum = (precipitation.where(cloud_types == 2) * time_interval).sum(dim='time')
        chunk_st_count_sum = (cloud_types == 2).sum(dim='time')
        
        # Non-deep convective (nd): cloud_types == 3
        chunk_nd_pcp_sum = (precipitation.where(cloud_types == 3) * time_interval).sum(dim='time')
        chunk_nd_count_sum = (cloud_types == 3).sum(dim='time')
        
        # Drizzle (dz): cloud_types == 4
        chunk_dz_pcp_sum = (precipitation.where(cloud_types == 4) * time_interval).sum(dim='time')
        chunk_dz_count_sum = (cloud_types == 4).sum(dim='time')
        
        # Accumulate results
        if totprecip_sum is None:
            # Initialize with first chunk results
            totprecip_sum = chunk_totprecip

            # Original features
            mcs_pcp_sum = chunk_mcs_pcp_sum
            mcs_count_sum = chunk_mcs_count_sum
            mcs_pcp_count_sum = chunk_mcs_pcp_count_sum

            ar_pcp_sum = chunk_ar_pcp_sum
            ar_count_sum = chunk_ar_count_sum
            ar_pcp_count_sum = chunk_ar_pcp_count_sum

            etc_pcp_sum = chunk_etc_pcp_sum
            etc_count_sum = chunk_etc_count_sum
            etc_pcp_count_sum = chunk_etc_pcp_count_sum

            tc_pcp_sum = chunk_tc_pcp_sum
            tc_count_sum = chunk_tc_count_sum
            tc_pcpcount_sum = chunk_tc_pcp_count_sum

            # Isolated features
            # MCS results
            mcs_iso_pcp_sum = chunk_mcs_iso_pcp_sum
            mcs_iso_count_sum = chunk_mcs_iso_count_sum
            mcs_iso_pcpcount_sum = chunk_mcs_iso_pcp_count_sum
            
            # AR results
            ar_iso_pcp_sum = chunk_ar_iso_pcp_sum
            ar_iso_count_sum = chunk_ar_iso_count_sum
            ar_iso_pcpcount_sum = chunk_ar_iso_pcp_count_sum
                       
            # ETC results
            etc_iso_pcp_sum = chunk_etc_iso_pcp_sum
            etc_iso_count_sum = chunk_etc_iso_count_sum
            etc_iso_pcpcount_sum = chunk_etc_iso_pcp_count_sum

            # Co-occurrence features
            mcs_ar_pcp_sum = chunk_mcs_ar_pcp_sum
            mcs_ar_count_sum = chunk_mcs_ar_count_sum
            mcs_ar_pcp_count_sum = chunk_mcs_ar_pcp_count_sum

            mcs_etc_pcp_sum = chunk_mcs_etc_pcp_sum
            mcs_etc_count_sum = chunk_mcs_etc_count_sum
            mcs_etc_pcp_count_sum = chunk_mcs_etc_pcp_count_sum

            ar_etc_pcp_sum = chunk_ar_etc_pcp_sum
            ar_etc_count_sum = chunk_ar_etc_count_sum
            ar_etc_pcp_count_sum = chunk_ar_etc_pcp_count_sum

            mcs_ar_etc_pcp_sum = chunk_mcs_ar_etc_pcp_sum
            mcs_ar_etc_count_sum = chunk_mcs_ar_etc_count_sum
            mcs_ar_etc_pcp_count_sum = chunk_mcs_ar_etc_pcp_count_sum
            
            # Cloud types
            dc_pcp_sum = chunk_dc_pcp_sum
            dc_count_sum = chunk_dc_count_sum
            
            st_pcp_sum = chunk_st_pcp_sum
            st_count_sum = chunk_st_count_sum
            
            nd_pcp_sum = chunk_nd_pcp_sum
            nd_count_sum = chunk_nd_count_sum
            
            dz_pcp_sum = chunk_dz_pcp_sum
            dz_count_sum = chunk_dz_count_sum
        else:
            # Add subsequent chunk results
            totprecip_sum += chunk_totprecip

            # Original features
            mcs_pcp_sum += chunk_mcs_pcp_sum
            mcs_count_sum += chunk_mcs_count_sum
            mcs_pcp_count_sum += chunk_mcs_pcp_count_sum

            ar_pcp_sum += chunk_ar_pcp_sum
            ar_count_sum += chunk_ar_count_sum
            ar_pcp_count_sum += chunk_ar_pcp_count_sum

            etc_pcp_sum += chunk_etc_pcp_sum
            etc_count_sum += chunk_etc_count_sum
            etc_pcp_count_sum += chunk_etc_pcp_count_sum

            tc_pcp_sum += chunk_tc_pcp_sum
            tc_count_sum += chunk_tc_count_sum
            tc_pcpcount_sum += chunk_tc_pcp_count_sum

            # Isolated features
            # MCS results
            mcs_iso_pcp_sum += chunk_mcs_iso_pcp_sum
            mcs_iso_count_sum += chunk_mcs_iso_count_sum
            mcs_iso_pcpcount_sum += chunk_mcs_iso_pcp_count_sum
            
            # AR results
            ar_iso_pcp_sum += chunk_ar_iso_pcp_sum
            ar_iso_count_sum += chunk_ar_iso_count_sum
            ar_iso_pcpcount_sum += chunk_ar_iso_pcp_count_sum            
            
            # ETC results
            etc_iso_pcp_sum += chunk_etc_iso_pcp_sum
            etc_iso_count_sum += chunk_etc_iso_count_sum
            etc_iso_pcpcount_sum += chunk_etc_iso_pcp_count_sum

            # Co-occurrence features
            mcs_ar_pcp_sum += chunk_mcs_ar_pcp_sum
            mcs_ar_count_sum += chunk_mcs_ar_count_sum
            mcs_ar_pcp_count_sum += chunk_mcs_ar_pcp_count_sum

            mcs_etc_pcp_sum += chunk_mcs_etc_pcp_sum
            mcs_etc_count_sum += chunk_mcs_etc_count_sum
            mcs_etc_pcp_count_sum += chunk_mcs_etc_pcp_count_sum

            ar_etc_pcp_sum += chunk_ar_etc_pcp_sum
            ar_etc_count_sum += chunk_ar_etc_count_sum
            ar_etc_pcp_count_sum += chunk_ar_etc_pcp_count_sum

            mcs_ar_etc_pcp_sum += chunk_mcs_ar_etc_pcp_sum
            mcs_ar_etc_count_sum += chunk_mcs_ar_etc_count_sum
            mcs_ar_etc_pcp_count_sum += chunk_mcs_ar_etc_pcp_count_sum

            # Cloud types
            dc_pcp_sum += chunk_dc_pcp_sum
            dc_count_sum += chunk_dc_count_sum
            
            st_pcp_sum += chunk_st_pcp_sum
            st_count_sum += chunk_st_count_sum
            
            nd_pcp_sum += chunk_nd_pcp_sum
            nd_count_sum += chunk_nd_count_sum
            
            dz_pcp_sum += chunk_dz_pcp_sum
            dz_count_sum += chunk_dz_count_sum

        # Explicitly delete chunk data to free memory
        del chunk_mcs_pcp_sum, chunk_mcs_count_sum, chunk_mcs_pcp_count_sum
        del chunk_ar_pcp_sum, chunk_ar_count_sum, chunk_ar_pcp_count_sum
        del chunk_etc_pcp_sum, chunk_etc_count_sum, chunk_etc_pcp_count_sum
        del chunk_tc_pcp_sum, chunk_tc_count_sum, chunk_tc_pcp_count_sum
        del chunk_ds, mcs_isolated_mask, ar_isolated_mask, etc_isolated_mask, tc_mask, precipitation
        del chunk_mcs_iso_pcp_sum, chunk_mcs_iso_count_sum, chunk_mcs_iso_pcp_count_sum
        del chunk_ar_iso_pcp_sum, chunk_ar_iso_count_sum, chunk_ar_iso_pcp_count_sum
        del chunk_etc_iso_pcp_sum, chunk_etc_iso_count_sum, chunk_etc_iso_pcp_count_sum
        del chunk_totprecip
        del chunk_mcs_ar_pcp_sum, chunk_mcs_ar_count_sum, chunk_mcs_ar_pcp_count_sum
        del chunk_mcs_etc_pcp_sum, chunk_mcs_etc_count_sum, chunk_mcs_etc_pcp_count_sum
        del chunk_ar_etc_pcp_sum, chunk_ar_etc_count_sum, chunk_ar_etc_pcp_count_sum
        del chunk_mcs_ar_etc_pcp_sum, chunk_mcs_ar_etc_count_sum, chunk_mcs_ar_etc_pcp_count_sum
        del mcs_ar_2way_mask, mcs_etc_2way_mask, ar_etc_2way_mask, mcs_ar_etc_3way_mask
        del chunk_dc_pcp_sum, chunk_dc_count_sum
        del chunk_st_pcp_sum, chunk_st_count_sum
        del chunk_nd_pcp_sum, chunk_nd_count_sum
        del chunk_dz_pcp_sum, chunk_dz_count_sum
        del cloud_types
        gc.collect()
    
    # Return all statistics in a dictionary
    return {
        'time': std_time,
        'ntimes': ntimes,
        'time_interval': time_interval,
        'totprecip': totprecip_sum,

        # Original features
        'mcs_precip': mcs_pcp_sum,
        'mcs_count': mcs_count_sum,
        'mcs_precip_count': mcs_pcp_count_sum,

        'ar_precip': ar_pcp_sum,
        'ar_count': ar_count_sum,
        'ar_precip_count': ar_pcp_count_sum,

        'etc_precip': etc_pcp_sum,
        'etc_count': etc_count_sum,
        'etc_precip_count': etc_pcp_count_sum,

        'tc_precip': tc_pcp_sum,
        'tc_count': tc_count_sum,
        'tc_precip_count': tc_pcpcount_sum,

        # Isolated features
        # MCS statistics
        'mcs_iso_precip': mcs_iso_pcp_sum,
        'mcs_iso_count': mcs_iso_count_sum,
        'mcs_iso_precip_count': mcs_iso_pcpcount_sum,
        
        # AR statistics
        'ar_iso_precip': ar_iso_pcp_sum,
        'ar_iso_count': ar_iso_count_sum,
        'ar_iso_precip_count': ar_iso_pcpcount_sum,
                
        # ETC statistics
        'etc_iso_precip': etc_iso_pcp_sum,
        'etc_iso_count': etc_iso_count_sum,
        'etc_iso_precip_count': etc_iso_pcpcount_sum,

        # Co-occurrence features (2-way)
        'mcs_ar_precip': mcs_ar_pcp_sum,
        'mcs_ar_count': mcs_ar_count_sum,
        'mcs_ar_precip_count': mcs_ar_pcp_count_sum,

        'mcs_etc_precip': mcs_etc_pcp_sum,
        'mcs_etc_count': mcs_etc_count_sum,
        'mcs_etc_precip_count': mcs_etc_pcp_count_sum,

        'ar_etc_precip': ar_etc_pcp_sum,
        'ar_etc_count': ar_etc_count_sum,
        'ar_etc_precip_count': ar_etc_pcp_count_sum,
        
        # Co-occurrence features (3-way)
        'mcs_ar_etc_precip': mcs_ar_etc_pcp_sum,
        'mcs_ar_etc_count': mcs_ar_etc_count_sum,
        'mcs_ar_etc_precip_count': mcs_ar_etc_pcp_count_sum,
        
        # Cloud types
        # Deep convection
        'dc_precip': dc_pcp_sum,
        'dc_count': dc_count_sum,
        
        # Stratiform
        'st_precip': st_pcp_sum,
        'st_count': st_count_sum,
        
        # Non-deep convective
        'nd_precip': nd_pcp_sum,
        'nd_count': nd_count_sum,
        
        # Drizzle
        'dz_precip': dz_pcp_sum,
        'dz_count': dz_count_sum,
    }


def write_netcdf(results, ds, output_filename, zoom, pcp_thresh, logger=None):
    """
    Write monthly precipitation statistics to a NetCDF file.
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    logger.info(f'Preparing data for output file: {output_filename}')
    
    # Prepare data for output dataset
    times = [r['time'] for r in results]
    ntimes_values = [r['ntimes'] for r in results]
    time_interval = results[0]['time_interval']
    
    # Extract total precipitation values
    totprecip_values = [r['totprecip'] for r in results]

    # Extract values for original features
    mcs_precip_values = [r['mcs_precip'] for r in results]
    mcs_count_values = [r['mcs_count'] for r in results]
    mcs_precip_count_values = [r['mcs_precip_count'] for r in results]

    ar_precip_values = [r['ar_precip'] for r in results]
    ar_count_values = [r['ar_count'] for r in results]
    ar_precip_count_values = [r['ar_precip_count'] for r in results]

    etc_precip_values = [r['etc_precip'] for r in results]
    etc_count_values = [r['etc_count'] for r in results]
    etc_precip_count_values = [r['etc_precip_count'] for r in results]

    # Extract values for isolated features
    mcs_iso_precip_values = [r['mcs_iso_precip'] for r in results]
    mcs_iso_count_values = [r['mcs_iso_count'] for r in results]
    mcs_iso_precip_count_values = [r['mcs_iso_precip_count'] for r in results]

    ar_iso_precip_values = [r['ar_iso_precip'] for r in results]
    ar_iso_count_values = [r['ar_iso_count'] for r in results]
    ar_iso_precip_count_values = [r['ar_iso_precip_count'] for r in results]

    etc_iso_precip_values = [r['etc_iso_precip'] for r in results]
    etc_iso_count_values = [r['etc_iso_count'] for r in results]
    etc_iso_precip_count_values = [r['etc_iso_precip_count'] for r in results]

    tcprecip_values = [r['tc_precip'] for r in results]
    tccount_values = [r['tc_count'] for r in results]
    tcpcpct_values = [r['tc_precip_count'] for r in results]

    # Extract values for co-occurrence features (2-way)
    mcs_ar_precip_values = [r['mcs_ar_precip'] for r in results]
    mcs_ar_count_values = [r['mcs_ar_count'] for r in results]
    mcs_ar_precip_count_values = [r['mcs_ar_precip_count'] for r in results]

    mcs_etc_precip_values = [r['mcs_etc_precip'] for r in results]
    mcs_etc_count_values = [r['mcs_etc_count'] for r in results]
    mcs_etc_precip_count_values = [r['mcs_etc_precip_count'] for r in results]

    ar_etc_precip_values = [r['ar_etc_precip'] for r in results]
    ar_etc_count_values = [r['ar_etc_count'] for r in results]
    ar_etc_precip_count_values = [r['ar_etc_precip_count'] for r in results]

    # Extract values for co-occurrence features (3-way)
    mcs_ar_etc_precip_values = [r['mcs_ar_etc_precip'] for r in results]
    mcs_ar_etc_count_values = [r['mcs_ar_etc_count'] for r in results]
    mcs_ar_etc_precip_count_values = [r['mcs_ar_etc_precip_count'] for r in results]

    # Extract values for cloud types
    # Deep convection
    dc_precip_values = [r['dc_precip'] for r in results]
    dc_count_values = [r['dc_count'] for r in results]
    
    # Stratiform
    st_precip_values = [r['st_precip'] for r in results]
    st_count_values = [r['st_count'] for r in results]
    
    # Non-deep convective
    nd_precip_values = [r['nd_precip'] for r in results]
    nd_count_values = [r['nd_count'] for r in results]
    
    # Drizzle
    dz_precip_values = [r['dz_precip'] for r in results]
    dz_count_values = [r['dz_count'] for r in results]

    # Create output variables
    var_dict = {
        'precipitation': (['time', 'cell'], np.stack([r.values for r in totprecip_values])),
        'ntimes': (['time'], np.array(ntimes_values)),

        # Original features
        'mcs_precipitation': (['time', 'cell'], np.stack([r.values for r in mcs_precip_values])),
        'mcs_count': (['time', 'cell'], np.stack([r.values for r in mcs_count_values])),
        'mcs_precipitation_count': (['time', 'cell'], np.stack([r.values for r in mcs_precip_count_values])),

        'ar_precipitation': (['time', 'cell'], np.stack([r.values for r in ar_precip_values])),
        'ar_count': (['time', 'cell'], np.stack([r.values for r in ar_count_values])),
        'ar_precipitation_count': (['time', 'cell'], np.stack([r.values for r in ar_precip_count_values])),

        'etc_precipitation': (['time', 'cell'], np.stack([r.values for r in etc_precip_values])),
        'etc_count': (['time', 'cell'], np.stack([r.values for r in etc_count_values])),
        'etc_precipitation_count': (['time', 'cell'], np.stack([r.values for r in etc_precip_count_values])),

        'tc_precipitation': (['time', 'cell'], np.stack([r.values for r in tcprecip_values])),
        'tc_count': (['time', 'cell'], np.stack([r.values for r in tccount_values])),
        'tc_precipitation_count': (['time', 'cell'], np.stack([r.values for r in tcpcpct_values])),
        
        # Isolated features
        'mcs_iso_precipitation': (['time', 'cell'], np.stack([r.values for r in mcs_iso_precip_values])),
        'mcs_iso_count': (['time', 'cell'], np.stack([r.values for r in mcs_iso_count_values])),
        'mcs_iso_precipitation_count': (['time', 'cell'], np.stack([r.values for r in mcs_iso_precip_count_values])),
        
        'ar_iso_precipitation': (['time', 'cell'], np.stack([r.values for r in ar_iso_precip_values])),
        'ar_iso_count': (['time', 'cell'], np.stack([r.values for r in ar_iso_count_values])),
        'ar_iso_precipitation_count': (['time', 'cell'], np.stack([r.values for r in ar_iso_precip_count_values])),
        
        'etc_iso_precipitation': (['time', 'cell'], np.stack([r.values for r in etc_iso_precip_values])),
        'etc_iso_count': (['time', 'cell'], np.stack([r.values for r in etc_iso_count_values])),
        'etc_iso_precipitation_count': (['time', 'cell'], np.stack([r.values for r in etc_iso_precip_count_values])),

        # Co-occurrence features (2-way)
        'mcs_ar_precipitation': (['time', 'cell'], np.stack([r.values for r in mcs_ar_precip_values])),
        'mcs_ar_count': (['time', 'cell'], np.stack([r.values for r in mcs_ar_count_values])),
        'mcs_ar_precipitation_count': (['time', 'cell'], np.stack([r.values for r in mcs_ar_precip_count_values])),

        'mcs_etc_precipitation': (['time', 'cell'], np.stack([r.values for r in mcs_etc_precip_values])),
        'mcs_etc_count': (['time', 'cell'], np.stack([r.values for r in mcs_etc_count_values])),
        'mcs_etc_precipitation_count': (['time', 'cell'], np.stack([r.values for r in mcs_etc_precip_count_values])),

        'ar_etc_precipitation': (['time', 'cell'], np.stack([r.values for r in ar_etc_precip_values])),
        'ar_etc_count': (['time', 'cell'], np.stack([r.values for r in ar_etc_count_values])),
        'ar_etc_precipitation_count': (['time', 'cell'], np.stack([r.values for r in ar_etc_precip_count_values])),

        # Co-occurrence features (3-way)
        'mcs_ar_etc_precipitation': (['time', 'cell'], np.stack([r.values for r in mcs_ar_etc_precip_values])),
        'mcs_ar_etc_count': (['time', 'cell'], np.stack([r.values for r in mcs_ar_etc_count_values])),
        'mcs_ar_etc_precipitation_count': (['time', 'cell'], np.stack([r.values for r in mcs_ar_etc_precip_count_values])),
        
        # Cloud types
        # Deep convection
        'dc_precipitation': (['time', 'cell'], np.stack([r.values for r in dc_precip_values])),
        'dc_count': (['time', 'cell'], np.stack([r.values for r in dc_count_values])),
        
        # Stratiform
        'st_precipitation': (['time', 'cell'], np.stack([r.values for r in st_precip_values])),
        'st_count': (['time', 'cell'], np.stack([r.values for r in st_count_values])),
        
        # Non-deep convective
        'nd_precipitation': (['time', 'cell'], np.stack([r.values for r in nd_precip_values])),
        'nd_count': (['time', 'cell'], np.stack([r.values for r in nd_count_values])),
        
        # Drizzle
        'dz_precipitation': (['time', 'cell'], np.stack([r.values for r in dz_precip_values])),
        'dz_count': (['time', 'cell'], np.stack([r.values for r in dz_count_values])),
    }
    
    # Create coordinates
    coord_dict = {
        # Use time objects directly, not their .values
        'time': times,
        'cell': (['cell'], ds['cell'].values),
        'lat': (['cell'], ds['lat'].values),
        'lon': (['cell'], ds['lon'].values),
    }
    
    # Add crs if available
    if 'crs' in ds:
        coord_dict['crs'] = ds['crs'].values
    
    # Create global attributes
    gattr_dict = {
        'Title': 'Monthly precipitation statistics by co-occurrence feature type',
        'contact': 'Zhe Feng, zhe.feng@pnnl.gov',
        # 'start_date': start_datetime,
        # 'end_date': end_datetime,
        'created_on': time.ctime(time.time()),
        'grid_type': 'HEALPix',
        'zoom_level': zoom,
        'time_interval': time_interval,
        'precipitation_threshold': pcp_thresh,
    }

    # Create output dataset
    dsout = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)

    # Add variable attributes
    dsout['cell'].attrs['long_name'] = 'HEALPix cell index'
    dsout['lon'].attrs['long_name'] = 'Longitude'
    dsout['lon'].attrs['units'] = 'degree'
    dsout['lat'].attrs['long_name'] = 'Latitude'
    dsout['lat'].attrs['units'] = 'degree'
    dsout['ntimes'].attrs['long_name'] = 'Number of hours during the month'
    dsout['ntimes'].attrs['units'] = 'count'
    
    # Total precipitation
    dsout['precipitation'].attrs['long_name'] = 'Total precipitation'
    dsout['precipitation'].attrs['units'] = 'mm'

    # Original features
    dsout['mcs_precipitation'].attrs['long_name'] = 'MCS precipitation (all)'
    dsout['mcs_precipitation'].attrs['units'] = 'mm'
    dsout['mcs_count'].attrs['long_name'] = 'Number of hours MCS is present'
    dsout['mcs_count'].attrs['units'] = 'hour'
    dsout['mcs_precipitation_count'].attrs['long_name'] = 'Number of hours MCS precipitation exceeds threshold'
    dsout['mcs_precipitation_count'].attrs['units'] = 'hour'

    dsout['ar_precipitation'].attrs['long_name'] = 'AR precipitation (all)'
    dsout['ar_precipitation'].attrs['units'] = 'mm'
    dsout['ar_count'].attrs['long_name'] = 'Number of hours AR is present'
    dsout['ar_count'].attrs['units'] = 'hour'
    dsout['ar_precipitation_count'].attrs['long_name'] = 'Number of hours AR precipitation exceeds threshold'
    dsout['ar_precipitation_count'].attrs['units'] = 'hour'

    dsout['etc_precipitation'].attrs['long_name'] = 'ETC precipitation (all)'
    dsout['etc_precipitation'].attrs['units'] = 'mm'
    dsout['etc_count'].attrs['long_name'] = 'Number of hours ETC is present'
    dsout['etc_count'].attrs['units'] = 'hour'
    dsout['etc_precipitation_count'].attrs['long_name'] = 'Number of hours ETC precipitation exceeds threshold'
    dsout['etc_precipitation_count'].attrs['units'] = 'hour'

    # TC attributes
    dsout['tc_precipitation'].attrs['long_name'] = 'TC precipitation'
    dsout['tc_precipitation'].attrs['units'] = 'mm'
    dsout['tc_count'].attrs['long_name'] = 'Number of hours TC is present'
    dsout['tc_count'].attrs['units'] = 'hour'
    dsout['tc_precipitation_count'].attrs['long_name'] = 'Number of hours TC precipitation exceeds threshold'
    dsout['tc_precipitation_count'].attrs['units'] = 'hour'
    
    # Isolated features
    # MCS attributes
    dsout['mcs_iso_precipitation'].attrs['long_name'] = 'MCS isolated precipitation'
    dsout['mcs_iso_precipitation'].attrs['units'] = 'mm'
    dsout['mcs_iso_count'].attrs['long_name'] = 'Number of hours MCS is present'
    dsout['mcs_iso_count'].attrs['units'] = 'hour'
    dsout['mcs_iso_precipitation_count'].attrs['long_name'] = 'Number of hours MCS isolated precipitation exceeds threshold'
    dsout['mcs_iso_precipitation_count'].attrs['units'] = 'hour'
    
    # AR attributes
    dsout['ar_iso_precipitation'].attrs['long_name'] = 'AR isolated precipitation'
    dsout['ar_iso_precipitation'].attrs['units'] = 'mm'
    dsout['ar_iso_count'].attrs['long_name'] = 'Number of hours AR is present'
    dsout['ar_iso_count'].attrs['units'] = 'hour'
    dsout['ar_iso_precipitation_count'].attrs['long_name'] = 'Number of hours AR isolated precipitation exceeds threshold'
    dsout['ar_iso_precipitation_count'].attrs['units'] = 'hour'
    
    # ETC attributes
    dsout['etc_iso_precipitation'].attrs['long_name'] = 'ETC isolated precipitation'
    dsout['etc_iso_precipitation'].attrs['units'] = 'mm'
    dsout['etc_iso_count'].attrs['long_name'] = 'Number of hours ETC is present'
    dsout['etc_iso_count'].attrs['units'] = 'hour'
    dsout['etc_iso_precipitation_count'].attrs['long_name'] = 'Number of hours ETC isolated precipitation exceeds threshold'
    dsout['etc_iso_precipitation_count'].attrs['units'] = 'hour'

    # Co-occurrence features (2-way)
    dsout['mcs_ar_precipitation'].attrs['long_name'] = 'MCS-AR co-occurrence precipitation'
    dsout['mcs_ar_precipitation'].attrs['units'] = 'mm'
    dsout['mcs_ar_count'].attrs['long_name'] = 'Number of hours MCS-AR co-occurrence is present'
    dsout['mcs_ar_count'].attrs['units'] = 'hour'
    dsout['mcs_ar_precipitation_count'].attrs['long_name'] = 'Number of hours MCS-AR co-occurrence precipitation exceeds threshold'
    dsout['mcs_ar_precipitation_count'].attrs['units'] = 'hour'

    dsout['mcs_etc_precipitation'].attrs['long_name'] = 'MCS-ETC co-occurrence precipitation'
    dsout['mcs_etc_precipitation'].attrs['units'] = 'mm'
    dsout['mcs_etc_count'].attrs['long_name'] = 'Number of hours MCS-ETC co-occurrence is present'
    dsout['mcs_etc_count'].attrs['units'] = 'hour'
    dsout['mcs_etc_precipitation_count'].attrs['long_name'] = 'Number of hours MCS-ETC co-occurrence precipitation exceeds threshold'
    dsout['mcs_etc_precipitation_count'].attrs['units'] = 'hour'

    dsout['ar_etc_precipitation'].attrs['long_name'] = 'AR-ETC co-occurrence precipitation'
    dsout['ar_etc_precipitation'].attrs['units'] = 'mm'
    dsout['ar_etc_count'].attrs['long_name'] = 'Number of hours AR-ETC co-occurrence is present'
    dsout['ar_etc_count'].attrs['units'] = 'hour'
    dsout['ar_etc_precipitation_count'].attrs['long_name'] = 'Number of hours AR-ETC co-occurrence precipitation exceeds threshold'
    dsout['ar_etc_precipitation_count'].attrs['units'] = 'hour'

    # Co-occurrence features (3-way)
    dsout['mcs_ar_etc_precipitation'].attrs['long_name'] = 'MCS-AR-ETC co-occurrence precipitation'
    dsout['mcs_ar_etc_precipitation'].attrs['units'] = 'mm'
    dsout['mcs_ar_etc_count'].attrs['long_name'] = 'Number of hours MCS-AR-ETC co-occurrence is present'
    dsout['mcs_ar_etc_count'].attrs['units'] = 'hour'
    dsout['mcs_ar_etc_precipitation_count'].attrs['long_name'] = 'Number of hours MCS-AR-ETC co-occurrence precipitation exceeds threshold'
    dsout['mcs_ar_etc_precipitation_count'].attrs['units'] = 'hour'

    # Cloud types
    # Deep convection
    dsout['dc_precipitation'].attrs['long_name'] = 'Deep convective cloud precipitation'
    dsout['dc_precipitation'].attrs['units'] = 'mm'
    dsout['dc_precipitation'].attrs['cloud_type_value'] = 1
    dsout['dc_count'].attrs['long_name'] = 'Number of hours deep convective cloud is present'
    dsout['dc_count'].attrs['units'] = 'hour'
    dsout['dc_count'].attrs['cloud_type_value'] = 1
    
    # Stratiform
    dsout['st_precipitation'].attrs['long_name'] = 'Stratiform cloud precipitation'
    dsout['st_precipitation'].attrs['units'] = 'mm'
    dsout['st_precipitation'].attrs['cloud_type_value'] = 2
    dsout['st_count'].attrs['long_name'] = 'Number of hours stratiform cloud is present'
    dsout['st_count'].attrs['units'] = 'hour'
    dsout['st_count'].attrs['cloud_type_value'] = 2
    
    # Non-deep convective
    dsout['nd_precipitation'].attrs['long_name'] = 'Non-deep convective cloud precipitation'
    dsout['nd_precipitation'].attrs['units'] = 'mm'
    dsout['nd_precipitation'].attrs['cloud_type_value'] = 3
    dsout['nd_count'].attrs['long_name'] = 'Number of hours non-deep convective cloud is present'
    dsout['nd_count'].attrs['units'] = 'hour'
    dsout['nd_count'].attrs['cloud_type_value'] = 3
    
    # Drizzle
    dsout['dz_precipitation'].attrs['long_name'] = 'Drizzle precipitation'
    dsout['dz_precipitation'].attrs['units'] = 'mm'
    dsout['dz_precipitation'].attrs['cloud_type_value'] = 4
    dsout['dz_count'].attrs['long_name'] = 'Number of hours drizzle cloud is present'
    dsout['dz_count'].attrs['units'] = 'hour'
    dsout['dz_count'].attrs['cloud_type_value'] = 4

    # Save the output file
    fillvalue = np.nan
    comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}

    logger.info(f'Writing output file: {output_filename}')
    dsout.to_netcdf(path=output_filename, mode='w', format='NETCDF4', 
                    unlimited_dims='time', encoding=encoding)
    logger.info(f'Successfully wrote: {output_filename}')
    
    return dsout


def subset_time_range(ds, start_datetime_str, end_datetime_str, logger=None):
    """
    Subset dataset to a time range, handling different calendar types robustly.
    
    Args:
        ds: xarray.Dataset
            Dataset to subset
        start_datetime_str: str
            Start datetime string in format 'YYYY-MM-DDTHH:MM'
        end_datetime_str: str
            End datetime string in format 'YYYY-MM-DDTHH:MM'
        logger: logging.Logger, optional
            Logger for status messages
            
    Returns:
        xarray.Dataset: Time-subsetted dataset
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    logger.info(f"Subsetting time range: {start_datetime_str} to {end_datetime_str}")
    
    # Parse datetime strings
    start_dt = pd.to_datetime(start_datetime_str)
    end_dt = pd.to_datetime(end_datetime_str)
    
    # Get the calendar type of the dataset
    time_values = ds.time.values
    calendar_type = type(time_values[0]).__name__
    
    logger.info(f"Dataset uses calendar type: {calendar_type}")
    
    if calendar_type in ['datetime64', 'Timestamp']:
        # Standard numpy datetime64 or pandas Timestamp
        logger.info("Using standard datetime subsetting")
        ds_subset = ds.sel(time=slice(start_datetime_str, end_datetime_str))
        
    elif 'cftime' in calendar_type.lower() or hasattr(time_values[0], 'calendar'):
        # cftime objects (e.g., DatetimeNoLeap, Datetime360Day, etc.)
        logger.info("Using cftime calendar subsetting")
        
        # Determine the specific cftime calendar
        sample_time = time_values[0]
        if hasattr(sample_time, 'calendar'):
            calendar_name = sample_time.calendar
        else:
            # Infer calendar from type name
            if 'NoLeap' in calendar_type:
                calendar_name = 'noleap'
            elif '360' in calendar_type:
                calendar_name = '360_day'
            elif 'Gregorian' in calendar_type:
                calendar_name = 'gregorian'
            else:
                calendar_name = 'standard'
        
        logger.info(f"Detected cftime calendar: {calendar_name}")
        
        # Convert start/end times to cftime objects with matching calendar
        try:
            start_cftime = cftime.datetime(
                start_dt.year, start_dt.month, start_dt.day,
                start_dt.hour, start_dt.minute, start_dt.second,
                calendar=calendar_name
            )
            end_cftime = cftime.datetime(
                end_dt.year, end_dt.month, end_dt.day,
                end_dt.hour, end_dt.minute, end_dt.second,
                calendar=calendar_name
            )
        except Exception as e:
            logger.warning(f"Failed to create cftime objects with calendar {calendar_name}: {e}")
            # Fallback: use the same type as the dataset
            start_cftime = type(sample_time)(
                start_dt.year, start_dt.month, start_dt.day,
                start_dt.hour, start_dt.minute, start_dt.second
            )
            end_cftime = type(sample_time)(
                end_dt.year, end_dt.month, end_dt.day,
                end_dt.hour, end_dt.minute, end_dt.second
            )
        
        # Create boolean mask for time selection
        time_mask = (ds.time >= start_cftime) & (ds.time <= end_cftime)
        ds_subset = ds.where(time_mask, drop=True)
        
    else:
        # Fallback: try converting dataset times to pandas datetime for comparison
        logger.warning(f"Unknown calendar type {calendar_type}, attempting fallback method")
        
        try:
            # Convert dataset times to pandas datetime for comparison
            if hasattr(time_values[0], 'year'):
                # Has year, month, day attributes (cftime-like)
                pd_times = [pd.Timestamp(t.year, t.month, t.day, 
                                       getattr(t, 'hour', 0), 
                                       getattr(t, 'minute', 0), 
                                       getattr(t, 'second', 0)) 
                           for t in time_values]
            else:
                # Try direct conversion
                pd_times = pd.to_datetime(time_values)
            
            # Create boolean mask
            time_mask = (pd.Series(pd_times) >= start_dt) & (pd.Series(pd_times) <= end_dt)
            ds_subset = ds.isel(time=time_mask)
            
        except Exception as e:
            logger.error(f"Fallback time subsetting failed: {e}")
            logger.error("Using string-based selection as last resort")
            ds_subset = ds.sel(time=slice(start_datetime_str, end_datetime_str))
    
    # Log results
    original_times = len(ds.time)
    subset_times = len(ds_subset.time)
    logger.info(f"Time subsetting complete: {original_times} -> {subset_times} time steps")
    
    if subset_times == 0:
        logger.warning("No time steps found in specified range!")
    else:
        first_time = ds_subset.time.values[0]
        last_time = ds_subset.time.values[-1]
        logger.info(f"Subset time range: {first_time} to {last_time}")
    
    return ds_subset

def main():
    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    # Start timing
    start_time = time.time()
    initial_memory = get_memory_usage()
    print(f"Initial memory usage: {initial_memory:.2f} GB")

    # Get the command-line arguments
    args_dict = parse_cmd_args()
    config_file = args_dict.get('config_file')
    catalog_source = args_dict.get('source')
    # zoom = args_dict.get('zoom')
    n_workers = args_dict.get('n_workers')
    threads_per_worker = args_dict.get('threads_per_worker')
    chunk_days = args_dict.get('chunk_days')
    pcp_thresh = args_dict.get('pcp_thresh')

    # chunk_days = 6
    # pcp_thresh = 0.1  # mm/h

    # Configuration parameters
    zoom = 8
    version = 'v1'
    parallel = True
    
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

    # Input combined mask file
    in_dir = "/pscratch/sd/w/wcmca1/hackathon/cof_masks/"
    in_basename = f"{source_name}_cofmasks_hp{zoom}_{version}.zarr"
    in_zarr = f"{in_dir}{in_basename}"

    output_dir = "/pscratch/sd/w/wcmca1/hackathon/cof_masks/stats/monthly/"
    os.makedirs(output_dir, exist_ok=True)
    output_filename = f"{output_dir}{source_name}_monthly_rainmap_cof_hp{zoom}_{version}.nc"

    logger.info(f"Using catalog source: {catalog_source}")
    logger.info(f"Source name: {source_name}")
    logger.info(f"Input file: {in_zarr}")
    logger.info(f"Output file: {output_filename}")


    # Setup Dask client
    client = setup_dask_client(parallel=parallel, n_workers=n_workers, threads_per_worker=threads_per_worker, logger=logger)

    try:
        # Load the mask dataset
        ds = xr.open_zarr(in_zarr, consolidated=True)
        ds = ds.pipe(egh.attach_coords)
        
        # Special treatment for certain datasets not in the catalog
        if catalog_source == "IR_IMERG":
            # Special case for IMERG data (not in catalog yet)
            dir_healpix = "/pscratch/sd/w/wcmca1/GPM/healpix/"
            in_basename = f"IMERG_V7_"
            time_res = "6H"
            in_zarr = f"{dir_healpix}{in_basename}{time_res}_zoom{zoom}_20190101_20211231.zarr"
            # Read IMERG dataset
            print(f"Loading IMERG dataset (NOT from catalog): {in_zarr}")
            ds_p = xr.open_zarr(in_zarr, consolidated=True)
            ds_p = ds_p.pipe(egh.attach_coords)

        elif catalog_source == "nicam_gl11":
            dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/nicam_gl11/shifted/"
            in_basename = f"NICAM_pr"
            time_res = "6h"
            in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
            # Read NICAM dataset
            print(f"Loading NICAM dataset (NOT from catalog): {in_zarr}")
            ds_p = xr.open_zarr(in_zarr, consolidated=True)
            ds_p = ds_p.pipe(egh.attach_coords)

        elif catalog_source == "um_glm_n2560_RAL3p3":
            dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/um_glm_n2560_RAL3p3/"
            in_basename = f"um_glm_n2560_RAL3p3_pr"
            time_res = "6h"
            in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
            # Read UM dataset
            print(f"Loading UM dataset (NOT from catalog): {in_zarr}")
            ds_p = xr.open_zarr(in_zarr, consolidated=True)
            ds_p = ds_p.pipe(egh.attach_coords)

        elif catalog_source == "casesm2_10km_nocumulus":
            dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/casesm2_10km_nocumulus/"
            in_basename = f"casesm2_10km_nocumulus_pr"
            time_res = "6h"
            in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
            # Read CASESM2 dataset
            print(f"Loading CASESM2 dataset (NOT from catalog): {in_zarr}")
            ds_p = xr.open_zarr(in_zarr, consolidated=True)
            ds_p = ds_p.pipe(egh.attach_coords)

        else:
            # Load the HEALPix catalog
            print(f"Loading HEALPix catalog: {catalog_file}")
            in_catalog = intake.open_catalog(catalog_file)
            if catalog_location:
                in_catalog = in_catalog[catalog_location]
            
            # Get the DataSet from the catalog
            ds_p = in_catalog[catalog_source](**catalog_params).to_dask()
            # Add lat/lon coordinates to the HEALPix DataSet
            ds_p = ds_p.pipe(egh.attach_coords)

        # Check liquid precipitaiton variable
        if varname_precip_liq in list(ds_p.keys()):
            # Convert liquid precipitation to mm/h
            pr = ds_p[varname_precip_liq] * pr_convert_factor
        # Check if the ice precipitation variable exist in the dataset
        if varname_precip_ice in list(ds_p.keys()):
            # Convert ice precipitation to liquid equivalent
            prs = ds_p[varname_precip_ice] * pr_convert_factor
            # Add ice precipitation to get total precipitation
            pr = pr + prs

        # Calendar conversion - check and convert calendars to match
        logger.info("Checking time coordinate calendars...")
        
        # Determine calendar types
        ds_p_calendar_type = type(ds_p.time.values[0]).__name__
        ds_calendar_type = type(ds.time.values[0]).__name__
        
        logger.info(f"Dataset calendars: ds_p uses {ds_p_calendar_type}, ds uses {ds_calendar_type}")
        
        # Convert ds time to match ds_p if they differ
        if ds_calendar_type != ds_p_calendar_type:
            logger.info(f"Converting ds time from {ds_calendar_type} to {ds_p_calendar_type}")
            
            # Get the calendar details from ds_p
            has_year_zero = True
            if hasattr(ds_p.time.values[0], 'has_year_zero'):
                has_year_zero = ds_p.time.values[0].has_year_zero
            
            # Convert datetime64 values to cftime DatetimeNoLeap objects
            new_times = []
            for t in ds.time.values:
                # Convert numpy datetime64 to pandas Timestamp to get date components
                pd_time = pd.Timestamp(t)
                # Create a matching cftime object
                dt_cftime = cftime.DatetimeNoLeap(
                    pd_time.year, pd_time.month, pd_time.day,
                    pd_time.hour, pd_time.minute, pd_time.second,
                    has_year_zero=has_year_zero
                )
                new_times.append(dt_cftime)
            
            # Create a new dataset with the converted time coordinate
            ds = ds.assign_coords(time=new_times)
            logger.info("Calendar conversion complete")

        # Find common time range across all datasets
        common_times = sorted(set(ds_p['time'].values)
                             .intersection(set(ds['time'].values)))
        if not common_times:
            logger.warning("No common time values between all datasets!")
            return None
        else:
            # Select only the common times in all datasets
            pr = pr.sel(time=common_times)
            ds = ds.sel(time=common_times)
            # Add precipitation to the dataset
            ds["pr"] = pr
            logger.info(f"Successfully merged datasets with {len(common_times)} common time points")

        # Subset to the specified time range using robust method
        if start_datetime and end_datetime:
            logger.info("Subsetting datasets to specified time range")
            ds = subset_time_range(ds, start_datetime, end_datetime, logger)
            
            if len(ds.time) == 0:
                logger.error("No data found in specified time range!")
                return None
        else:
            logger.info("No time range specified, using all available data")

        # Group by month and apply the processing function
        monthly_groups = ds.resample(time='1MS')

        # Check if client exists for parallel processing
        if client is not None:
            # Parallel processing with Dask
            logger.info("Running in parallel mode with Dask")
            delayed_results = []
            for month_start, month_ds in monthly_groups:
                print(f"Processing month: {month_start}")
                # print(f"Processing month: {month_start.strftime('%Y-%m')}")
                # Submit the processing job to the dask cluster
                delayed_result = client.submit(process_month_chunked, month_ds, chunk_days=chunk_days, pcp_thresh=pcp_thresh)
                delayed_results.append(delayed_result)
            
            # Clear line and show progress tracking
            print("\nTracking progress of all months processing in parallel:")
            from dask.distributed import progress
            progress(delayed_results)
            
            # Gather results (will wait for completion)
            results = client.gather(delayed_results)
        else:
            # Serial processing
            logger.info("Running in serial mode")
            results = []
            for month_start, month_ds in monthly_groups:
                print(f"Processing month: {month_start}")
                # print(f"Processing month: {month_start.strftime('%Y-%m')}")
                # Process directly without Dask
                result = process_month_chunked(month_ds, chunk_days=chunk_days, pcp_thresh=pcp_thresh)
                results.append(result)

        # Write output to NetCDF file
        write_netcdf(results, ds, output_filename, zoom, pcp_thresh, logger)

        # Calculate total processing time
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        # Store success info to print after Dask cleanup
        success_info = {
            'elapsed_time': elapsed_time,
            'output_filename': output_filename,
            'n_months': len(results),
        }

    finally:
        # Always cleanup client
        if client and parallel:
            # Suppress Dask shutdown messages by temporarily raising log level
            logging.getLogger('distributed').setLevel(logging.CRITICAL)
            logging.getLogger('distributed.worker').setLevel(logging.CRITICAL)
            logging.getLogger('distributed.nanny').setLevel(logging.CRITICAL)
            
            logger.info("Shutting down Dask client")
            client.close()
        
        # Print success message after Dask cleanup (so it's always visible at the end)
        if 'success_info' in locals():
            elapsed = success_info['elapsed_time']
            print(f"\n{'='*80}")
            print(f"✅ PROCESSING COMPLETE!")
            print(f"{'='*80}")
            print(f"Output: {success_info['output_filename']}")
            print(f"Months processed: {success_info['n_months']}")
            print(f"Total time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
            print(f"{'='*80}\n")


if __name__ == "__main__":
    main()