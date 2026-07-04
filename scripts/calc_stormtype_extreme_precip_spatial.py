#!/usr/bin/env python
"""
Calculate spatial distribution of extreme precipitation by storm type using Dask parallelization.

This script processes precipitation and storm mask data to attribute extreme precipitation
events to specific storm types (MCS, AR, ETC, TC, co-occurrences, and cloud types).
The processing is parallelized over the time dimension using Dask for efficiency.

Author: Zhe Feng (zhe.feng@pnnl.gov)
Date: December 2025
"""

import os
import sys
import argparse
import time
import numpy as np
import xarray as xr
import pandas as pd
import cftime
import yaml
from pathlib import Path
import warnings
import logging
import easygems.healpix as egh
import intake
from dask.distributed import Client, progress
import dask
from dask.diagnostics import ProgressBar

warnings.filterwarnings('ignore')


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


def load_precipitation_data(catalog_file, catalog_source, config):
    """
    Load precipitation data from catalog or direct zarr path.
    
    Args:
        catalog_file: str
            URL to the catalog file
        catalog_source: str
            Source name in the catalog
        config: dict
            Configuration dictionary
    
    Returns:
        xr.DataArray: Precipitation data (mm/h)
    """
    catalog_location = config.get('catalog_location', 'NERSC')
    catalog_params = config.get('catalog_params', {}).copy()
    varname_precip_liq = config.get('varname_precip_liq')
    varname_precip_ice = config.get('varname_precip_ice')
    pr_convert_factor = config.get('pr_convert_factor')
    zoom = 8
    
    # Special treatment for certain datasets not in catalog
    if catalog_source == "IR_IMERG":
        dir_healpix = "/pscratch/sd/w/wcmca1/GPM/healpix/"
        in_basename = f"IMERG_V7_"
        time_res = "6H"
        in_zarr = f"{dir_healpix}{in_basename}{time_res}_zoom{zoom}_20190101_20211231.zarr"
        print(f"Loading IMERG dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)
    
    elif catalog_source == "scream_ne120":
        dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/scream/"
        in_basename = f"scream_pr"
        time_res = "6h"
        in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
        print(f"Loading SCREAM dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)
    
    elif catalog_source == "nicam_gl11":
        dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/nicam_gl11/shifted/"
        in_basename = f"NICAM_pr"
        time_res = "6h"
        in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
        print(f"Loading NICAM dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)
    
    elif catalog_source == "um_glm_n2560_RAL3p3":
        dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/um_glm_n2560_RAL3p3/"
        in_basename = f"um_glm_n2560_RAL3p3_pr"
        time_res = "6h"
        in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
        print(f"Loading UM dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)
    
    elif catalog_source == "casesm2_10km_nocumulus":
        dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/casesm2_10km_nocumulus/"
        in_basename = f"casesm2_10km_nocumulus_pr"
        time_res = "6h"
        in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
        print(f"Loading CASESM2 dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)
    
    else:
        # Load from catalog
        print(f"Loading HEALPix catalog: {catalog_file}")
        in_catalog = intake.open_catalog(catalog_file)
        if catalog_location:
            in_catalog = in_catalog[catalog_location]
        
        ds_p = in_catalog[catalog_source](**catalog_params).to_dask()
        ds_p = ds_p.pipe(egh.attach_coords)
    
    # Get precipitation variable
    if varname_precip_liq in list(ds_p.keys()):
        pr = ds_p[varname_precip_liq] * pr_convert_factor
    
    # Add ice precipitation if available
    if varname_precip_ice in list(ds_p.keys()):
        prs = ds_p[varname_precip_ice] * pr_convert_factor
        pr = pr + prs
    
    return pr


def create_union_mask(mask_list):
    """
    Create a union mask from multiple mask variables.
    
    Parameters:
    -----------
    mask_list : list of xarray.DataArray
        List of mask DataArrays to combine
    
    Returns:
    --------
    xarray.DataArray
        Binary mask representing the union of all input masks
    """
    if not mask_list:
        raise ValueError("mask_list cannot be empty")
    
    binary_masks = [(mask > 0).astype(float) for mask in mask_list]
    summed_mask = sum(binary_masks)
    union_mask = (summed_mask > 0).astype(float)
    
    return union_mask


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", 
        level=logging.INFO
    )


def standardize_calendar_type(ds, reference_time_values, logger=None):
    """
    Convert dataset time coordinate to match a reference calendar type.
    
    This function handles conversion between different calendar types (e.g., cftime, 
    datetime64, pandas Timestamp) to ensure compatibility between datasets.
    
    Parameters:
    -----------
    ds : xr.Dataset or xr.DataArray
        Dataset whose time coordinate needs conversion
    reference_time_values : np.ndarray
        Reference time values whose calendar type should be matched
    logger : logging.Logger, optional
        Logger instance for status messages
    
    Returns:
    --------
    xr.Dataset or xr.DataArray : Dataset with converted time coordinate
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Determine calendar types
    ref_calendar_type = type(reference_time_values[0]).__name__
    ds_calendar_type = type(ds.time.values[0]).__name__
    
    logger.info(f"Reference calendar: {ref_calendar_type}")
    logger.info(f"Dataset calendar: {ds_calendar_type}")
    
    # Convert ds time to match reference if they differ
    if ds_calendar_type != ref_calendar_type:
        logger.info(f"⚠️  Calendar mismatch detected! Converting dataset calendar to match reference...")
        logger.info(f"   Converting from {ds_calendar_type} to {ref_calendar_type}")
        
        # Get the calendar details from reference
        has_year_zero = True
        if hasattr(reference_time_values[0], 'has_year_zero'):
            has_year_zero = reference_time_values[0].has_year_zero
        
        # Convert datetime64 values to cftime objects or vice versa
        new_times = []
        for t in ds.time.values:
            if 'datetime64' in ref_calendar_type or ref_calendar_type == 'Timestamp':
                # Convert to datetime64/Timestamp
                if hasattr(t, 'year'):
                    # From cftime to datetime64
                    pd_time = pd.Timestamp(t.year, t.month, t.day, 
                                         getattr(t, 'hour', 0), 
                                         getattr(t, 'minute', 0), 
                                         getattr(t, 'second', 0))
                else:
                    # Already datetime64
                    pd_time = pd.Timestamp(t)
                new_times.append(pd_time)
            else:
                # Convert to cftime (same type as reference)
                if hasattr(t, 'year'):
                    # Already cftime, convert to pandas then to target cftime
                    pd_time = pd.Timestamp(t.year, t.month, t.day,
                                         getattr(t, 'hour', 0),
                                         getattr(t, 'minute', 0), 
                                         getattr(t, 'second', 0))
                else:
                    # From datetime64 to cftime
                    pd_time = pd.Timestamp(t)
                
                # Create matching cftime object
                if 'NoLeap' in ref_calendar_type or ref_calendar_type == 'DatetimeNoLeap':
                    dt_cftime = cftime.DatetimeNoLeap(
                        pd_time.year, pd_time.month, pd_time.day,
                        pd_time.hour, pd_time.minute, pd_time.second,
                        has_year_zero=has_year_zero
                    )
                elif '360' in ref_calendar_type or ref_calendar_type == 'Datetime360Day':
                    dt_cftime = cftime.Datetime360Day(
                        pd_time.year, pd_time.month, pd_time.day,
                        pd_time.hour, pd_time.minute, pd_time.second,
                        has_year_zero=has_year_zero
                    )
                else:
                    # Default to standard calendar
                    dt_cftime = cftime.datetime(
                        pd_time.year, pd_time.month, pd_time.day,
                        pd_time.hour, pd_time.minute, pd_time.second
                    )
                new_times.append(dt_cftime)
        
        # Create a new dataset with the converted time coordinate
        ds = ds.assign_coords(time=new_times)
        logger.info("✅ Calendar conversion complete")
    else:
        logger.info("✅ Calendars match, no conversion needed")
    
    return ds


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


def process_single_timestep(pr_t, ds_t, pr_threshold, compute_cloud_types=True):
    """
    Process a single time step to count extreme precipitation by storm type.
    This function is designed to be called in parallel via Dask.
    
    Parameters:
    -----------
    pr_t : xr.DataArray
        Precipitation for single time step (cell,)
    ds_t : xr.Dataset
        Storm masks for single time step
    pr_threshold : float or xr.DataArray
        Precipitation threshold
    compute_cloud_types : bool
        Whether to process cloud types
    
    Returns:
    --------
    dict : Dictionary with spatial arrays for each storm type (counts and precip amounts)
    """
    # Create union masks
    mcs_ar_2way = create_union_mask([
        ds_t['mcs_ar_overlap_mask'],
        ds_t['ar_mcs_overlap_mask'],
    ])
    mcs_etc_2way = create_union_mask([
        ds_t['mcs_etc_overlap_mask'],
        ds_t['etc_mcs_overlap_mask'],
    ])
    ar_etc_2way = create_union_mask([
        ds_t['ar_etc_overlap_mask'],
        ds_t['etc_ar_overlap_mask'],
    ])
    mcs_ar_etc_3way = create_union_mask([
        ds_t['mcs_ar_etc_overlap_mask'],
        ds_t['ar_mcs_etc_overlap_mask'],
        ds_t['etc_mcs_ar_overlap_mask'],
    ])
    
    # Identify extreme precipitation cells
    extreme_mask = pr_t > pr_threshold
    
    # Initialize results
    results = {
        'extreme_mask': extreme_mask.astype(float),
    }
    
    # Track assigned cells
    assigned_mask = xr.zeros_like(extreme_mask, dtype=bool)
    
    # Priority order: 3-way > 2-way > isolated > TC > cloud types
    
    # 1. Three-way co-occurrence
    mask_3way = mcs_ar_etc_3way > 0
    extreme_3way_mask = (extreme_mask & mask_3way & ~assigned_mask)
    results['mcs_ar_etc_3way'] = extreme_3way_mask.astype(float)
    results['mcs_ar_etc_3way_pr'] = xr.where(extreme_3way_mask, pr_t, 0.0)
    assigned_mask = assigned_mask | extreme_3way_mask
    
    # 2. Two-way co-occurrences
    mask_mcs_ar = mcs_ar_2way > 0
    extreme_mcs_ar_mask = (extreme_mask & mask_mcs_ar & ~assigned_mask)
    results['mcs_ar_2way'] = extreme_mcs_ar_mask.astype(float)
    results['mcs_ar_2way_pr'] = xr.where(extreme_mcs_ar_mask, pr_t, 0.0)
    assigned_mask = assigned_mask | extreme_mcs_ar_mask
    
    mask_mcs_etc = mcs_etc_2way > 0
    extreme_mcs_etc_mask = (extreme_mask & mask_mcs_etc & ~assigned_mask)
    results['mcs_etc_2way'] = extreme_mcs_etc_mask.astype(float)
    results['mcs_etc_2way_pr'] = xr.where(extreme_mcs_etc_mask, pr_t, 0.0)
    assigned_mask = assigned_mask | extreme_mcs_etc_mask
    
    mask_ar_etc = ar_etc_2way > 0
    extreme_ar_etc_mask = (extreme_mask & mask_ar_etc & ~assigned_mask)
    results['ar_etc_2way'] = extreme_ar_etc_mask.astype(float)
    results['ar_etc_2way_pr'] = xr.where(extreme_ar_etc_mask, pr_t, 0.0)
    assigned_mask = assigned_mask | extreme_ar_etc_mask
    
    # 3. Isolated storms
    mask_mcs_iso = ds_t['mcs_isolated_mask'] > 0
    extreme_mcs_iso_mask = (extreme_mask & mask_mcs_iso & ~assigned_mask)
    results['mcs_isolated'] = extreme_mcs_iso_mask.astype(float)
    results['mcs_isolated_pr'] = xr.where(extreme_mcs_iso_mask, pr_t, 0.0)
    assigned_mask = assigned_mask | extreme_mcs_iso_mask
    
    mask_ar_iso = ds_t['ar_isolated_mask'] > 0
    extreme_ar_iso_mask = (extreme_mask & mask_ar_iso & ~assigned_mask)
    results['ar_isolated'] = extreme_ar_iso_mask.astype(float)
    results['ar_isolated_pr'] = xr.where(extreme_ar_iso_mask, pr_t, 0.0)
    assigned_mask = assigned_mask | extreme_ar_iso_mask
    
    mask_etc_iso = ds_t['etc_isolated_mask'] > 0
    extreme_etc_iso_mask = (extreme_mask & mask_etc_iso & ~assigned_mask)
    results['etc_isolated'] = extreme_etc_iso_mask.astype(float)
    results['etc_isolated_pr'] = xr.where(extreme_etc_iso_mask, pr_t, 0.0)
    assigned_mask = assigned_mask | extreme_etc_iso_mask
    
    # 4. TC
    mask_tc = ds_t['tc_mask'] > 0
    extreme_tc_mask = (extreme_mask & mask_tc & ~assigned_mask)
    results['tc'] = extreme_tc_mask.astype(float)
    results['tc_pr'] = xr.where(extreme_tc_mask, pr_t, 0.0)
    assigned_mask = assigned_mask | extreme_tc_mask
    
    # 5. Cloud types
    if compute_cloud_types and 'cloud_types' in ds_t:
        cloud_types = ds_t['cloud_types']
        
        # Remove tracked storms from cloud types
        cloud_types_cleaned = xr.where(
            ((ds_t['mcs_isolated_mask'] > 0) |
             (ds_t['ar_isolated_mask'] > 0) |
             (ds_t['etc_isolated_mask'] > 0) |
             (ds_t['tc_mask'] > 0) |
             (mcs_ar_2way > 0) |
             (mcs_etc_2way > 0) |
             (ar_etc_2way > 0) |
             (mcs_ar_etc_3way > 0)), 
            0, cloud_types, keep_attrs=False
        )
        
        # DC = 1
        mask_dc = (cloud_types_cleaned == 1)
        extreme_dc_mask = (extreme_mask & mask_dc & ~assigned_mask)
        results['dc'] = extreme_dc_mask.astype(float)
        results['dc_pr'] = xr.where(extreme_dc_mask, pr_t, 0.0)
        assigned_mask = assigned_mask | extreme_dc_mask
        
        # ND = 3 (non-deep convective; matches upstream cloud-type convention in
        # make_mcs_swath_masks.py: 1=deep convective, 2=stratiform, 3=non-deep convective, 4=drizzle)
        mask_nd = (cloud_types_cleaned == 3)
        extreme_nd_mask = (extreme_mask & mask_nd & ~assigned_mask)
        results['nd'] = extreme_nd_mask.astype(float)
        results['nd_pr'] = xr.where(extreme_nd_mask, pr_t, 0.0)
        assigned_mask = assigned_mask | extreme_nd_mask

        # ST = 2 (stratiform)
        mask_st = (cloud_types_cleaned == 2)
        extreme_st_mask = (extreme_mask & mask_st & ~assigned_mask)
        results['st'] = extreme_st_mask.astype(float)
        results['st_pr'] = xr.where(extreme_st_mask, pr_t, 0.0)
        assigned_mask = assigned_mask | extreme_st_mask
        
        # DZ = 4
        mask_dz = (cloud_types_cleaned == 4)
        extreme_dz_mask = (extreme_mask & mask_dz & ~assigned_mask)
        results['dz'] = extreme_dz_mask.astype(float)
        results['dz_pr'] = xr.where(extreme_dz_mask, pr_t, 0.0)
        assigned_mask = assigned_mask | extreme_dz_mask
    else:
        results['dc'] = xr.zeros_like(extreme_mask, dtype=float)
        results['nd'] = xr.zeros_like(extreme_mask, dtype=float)
        results['st'] = xr.zeros_like(extreme_mask, dtype=float)
        results['dz'] = xr.zeros_like(extreme_mask, dtype=float)
        results['dc_pr'] = xr.zeros_like(pr_t)
        results['nd_pr'] = xr.zeros_like(pr_t)
        results['st_pr'] = xr.zeros_like(pr_t)
        results['dz_pr'] = xr.zeros_like(pr_t)
    
    # Unassigned
    unassigned_mask = (extreme_mask & ~assigned_mask)
    results['unassigned'] = unassigned_mask.astype(float)
    results['unassigned_pr'] = xr.where(unassigned_mask, pr_t, 0.0)
    
    # Total extreme precipitation amount
    results['total_extreme_pr'] = xr.where(extreme_mask, pr_t, 0.0)
    
    return results


def process_timeseries_dask(pr, pr_threshold, ds, percentile_name='P90',
                             compute_cloud_types=True, n_workers=8, batch_size=200):
    """
    Process all time steps using Dask for parallel processing.
    
    Parameters:
    -----------
    pr : xr.DataArray
        Precipitation data (time, cell)
    pr_threshold : float or xr.DataArray
        Precipitation threshold
    ds : xr.Dataset
        Storm masks dataset (time, cell)
    percentile_name : str
        Name of percentile (e.g., 'P90')
    compute_cloud_types : bool
        Whether to compute cloud type contributions
    n_workers : int
        Number of Dask workers
    batch_size : int
        Number of time steps to process per batch. Reduce if running out of memory.
    
    Returns:
    --------
    xr.Dataset : Spatial counts dataset
    """
    n_times = len(pr.time)
    n_batches = (n_times + batch_size - 1) // batch_size
    print(f"\nProcessing {n_times} time steps with Dask ({n_workers} workers)...")
    print(f"Percentile threshold: {percentile_name}")
    print(f"Batch size: {batch_size} time steps ({n_batches} batches total)")
    
    # Initialize accumulators (one set of arrays for the whole run)
    template = pr.isel(time=0).compute()
    storm_type_keys = ['mcs_isolated', 'ar_isolated', 'etc_isolated', 'tc',
                       'mcs_ar_2way', 'mcs_etc_2way', 'ar_etc_2way', 'mcs_ar_etc_3way',
                       'dc', 'nd', 'st', 'dz', 'unassigned']
    
    accumulated_counts = {k: xr.zeros_like(template) for k in ['total_extreme'] + storm_type_keys}
    accumulated_precip = {k: xr.zeros_like(template) for k in ['total_extreme'] + storm_type_keys}
    
    # Process in batches to avoid holding all results in memory simultaneously
    for batch_idx in range(n_batches):
        t_start = batch_idx * batch_size
        t_end = min(t_start + batch_size, n_times)
        print(f"\nBatch {batch_idx + 1}/{n_batches}: time steps {t_start}–{t_end - 1}")
        
        delayed_results = []
        for t in range(t_start, t_end):
            pr_t = pr.isel(time=t)
            ds_t = ds.isel(time=t)
            delayed_result = dask.delayed(process_single_timestep)(
                pr_t, ds_t, pr_threshold, compute_cloud_types
            )
            delayed_results.append(delayed_result)
        
        with ProgressBar():
            batch_results = dask.compute(*delayed_results, scheduler='threads', num_workers=n_workers)
        
        # Accumulate and immediately discard batch results
        for result in batch_results:
            accumulated_counts['total_extreme'] += result['extreme_mask']
            accumulated_precip['total_extreme'] += result['total_extreme_pr']
            for key in storm_type_keys:
                accumulated_counts[key] += result[key]
                accumulated_precip[key] += result[f'{key}_pr']
        del batch_results
    
    print("\n✅ All batches complete!")
    print("Calculating precipitation fractions...")
    
    # Create output dataset
    storm_type_order = [
        'mcs_isolated', 'ar_isolated', 'etc_isolated', 'tc',
        'mcs_ar_2way', 'mcs_etc_2way', 'ar_etc_2way', 'mcs_ar_etc_3way',
        'dc', 'nd', 'st', 'dz', 'unassigned'
    ]
    storm_type_names = [
        'MCS (isolated)', 'AR (isolated)', 'ETC (isolated)', 'TC',
        'MCS+AR', 'MCS+ETC', 'AR+ETC', 'MCS+AR+ETC',
        'Deep Convection', 'Non-deep Convection', 'Stratiform', 'Drizzle',
        'Unassigned'
    ]
    
    # Calculate precipitation fractions
    # Fraction = sum(stormtype_extreme_precip) / sum(all_extreme_precip)
    accumulated_fractions = {}
    total_extreme_precip = accumulated_precip['total_extreme']
    
    for key in storm_type_order:
        # Calculate fraction: avoid division by zero
        frac = xr.where(
            total_extreme_precip > 0,
            accumulated_precip[key] / total_extreme_precip,
            0.0
        )
        accumulated_fractions[key] = frac
    
    print("✅ Fraction calculation complete!")
    
    # Create dataset
    ds_out = xr.Dataset(coords={'cell': accumulated_counts['total_extreme'].cell})
    
    # Add total extreme count and precipitation
    ds_out['total_extreme_count'] = accumulated_counts['total_extreme']
    ds_out['total_extreme_count'].attrs = {
        'long_name': 'Total count of extreme precipitation occurrences',
        'units': 'count',
        'percentile': percentile_name
    }
    
    ds_out['total_extreme_precip'] = accumulated_precip['total_extreme']
    ds_out['total_extreme_precip'].attrs = {
        'long_name': 'Total extreme precipitation amount',
        'units': 'mm/h',
        'percentile': percentile_name,
        'description': (
            'Sum of instantaneous precipitation rate (mm/h) over extreme time steps '
            '(pr > threshold). This is a sum over samples, not time-integrated by the '
            'sampling interval, so it is not a physical accumulated depth.'
        )
    }
    
    # Add counts and fractions for each storm type
    for i, (key, name) in enumerate(zip(storm_type_order, storm_type_names)):
        # Counts
        ds_out[f'{key}_count'] = accumulated_counts[key]
        ds_out[f'{key}_count'].attrs = {
            'long_name': f'Count of extreme precip from {name}',
            'units': 'count',
            'storm_type_code': i + 1
        }
        
        # Fractions
        ds_out[f'{key}_frac'] = accumulated_fractions[key]
        ds_out[f'{key}_frac'].attrs = {
            'long_name': f'Precipitation fraction from {name}',
            'units': 'fraction',
            'storm_type_code': i + 1,
            'description': 'sum(stormtype_extreme_precip) / sum(all_extreme_precip)'
        }
        # Add lat/lon from mask dataset (ensures consistent coordinates across sources)
    # Note: Using ds (mask) rather than pr because mask coordinates are standardized
    # and consistent across all sources, avoiding the 2D (source, cell) issue when
    # combining multiple sources with xr.open_mfdataset
    if 'lat' in ds.coords:
        ds_out['lat'] = ds.lat
        ds_out['lon'] = ds.lon
    
    # Global attributes
    ds_out.attrs = {
        'title': f'Spatial distribution of extreme precipitation by storm type ({percentile_name})',
        'percentile': percentile_name,
        'description': 'Per-cell counts and fractions of extreme precipitation by storm type. '
                      'Fractions calculated as: sum(stormtype_extreme_precip) / sum(all_extreme_precip)',
        'ntimes_processed': len(pr.time),
    }
    
    return ds_out


def save_spatial_results(ds_spatial, output_file, source_name, 
                          start_date, end_date, percentile_name):
    """
    Save spatial results to NetCDF with compression.
    
    Parameters:
    -----------
    ds_spatial : xr.Dataset
        Spatial results dataset
    output_file : str
        Output filename
    source_name : str
        Data source name
    start_date : str
        Start date
    end_date : str
        End date
    percentile_name : str
        Percentile name
    """
    # Update attributes
    ds_spatial.attrs.update({
        'title': f'Spatial storm type attribution for {percentile_name} extreme precipitation',
        'source_name': source_name,
        'start_date': start_date,
        'end_date': end_date,
        'percentile': percentile_name,
        'created_on': time.ctime(time.time()),
        'contact': 'Zhe Feng, zhe.feng@pnnl.gov',
        'description': 'Per-cell counts and fractions of extreme precipitation by storm type.'
    })

    # Setup encoding
    encoding = {}
    for var in ds_spatial.data_vars:
        if 'count' in var or var == 'total_extreme':
            encoding[var] = {'dtype': 'int32', 'zlib': True, 'complevel': 4}
    
    # Save
    print(f"\nSaving results to: {output_file}")
    ds_spatial.to_netcdf(output_file, encoding=encoding)
    print(f"✅ Saved!")
    
    # Print summary
    print(f"\n{'='*60}")
    print("Summary Statistics:")
    print(f"{'='*60}")
    print(f"Total cells with extreme precipitation: {int((ds_spatial['total_extreme_count'] > 0).sum().values):,}")
    print(f"Mean extreme precip count per cell: {float(ds_spatial['total_extreme_count'].mean().values):.2f}")
    print(f"Max extreme precip count: {int(ds_spatial['total_extreme_count'].max().values)}")
    print(f"Total extreme precipitation amount: {float(ds_spatial['total_extreme_precip'].sum().values):.2f} mm/h")
    
    # Storm type fraction statistics (global averages)
    print(f"\nGlobal Precipitation Fractions by Storm Type:")
    storm_type_labels = [
        'MCS (iso)', 'AR (iso)', 'ETC (iso)', 'TC',
        'MCS+AR', 'MCS+ETC', 'AR+ETC', 'MCS+AR+ETC',
        'DC', 'ND', 'ST', 'DZ', 'Unassigned'
    ]
    storm_type_keys = [
        'mcs_isolated', 'ar_isolated', 'etc_isolated', 'tc',
        'mcs_ar_2way', 'mcs_etc_2way', 'ar_etc_2way', 'mcs_ar_etc_3way',
        'dc', 'nd', 'st', 'dz', 'unassigned'
    ]
    
    total_extreme_precip = float(ds_spatial['total_extreme_precip'].sum().values)
    for label, key in zip(storm_type_labels, storm_type_keys):
        # Calculate global fraction from precipitation amounts
        precip_key = f'{key}_frac'
        if precip_key in ds_spatial:
            # Global average (weighted by total extreme precip at each cell)
            mean_frac = float(ds_spatial[precip_key].mean().values)
            print(f"  {label:20s}: {mean_frac*100:5.2f}%")
    print(f"{'='*60}\n")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Calculate spatial distribution of extreme precipitation by storm type',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--catalog_source', type=str, required=True,
                       help='Catalog source name (e.g., scream_ne120, IR_IMERG)')
    
    parser.add_argument('--config_file', type=str,
                       default='/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources.yaml',
                       help='Path to configuration YAML file')
    
    parser.add_argument('--percentiles', type=str, nargs='+',
                       default=['P90'],
                       help='Percentile names to process (e.g., P90 P95 P99)')
    
    parser.add_argument('--start_date', type=str, default=None,
                       help='Start date (YYYY-MM-DD). If not specified, uses full dataset')
    
    parser.add_argument('--end_date', type=str, default=None,
                       help='End date (YYYY-MM-DD). If not specified, uses full dataset')
    
    parser.add_argument('--output_dir', type=str,
                       default='/pscratch/sd/w/wcmca1/hackathon/extreme_precip',
                       help='Output directory for NetCDF files')
    
    parser.add_argument('--n_workers', type=int, default=8,
                       help='Number of Dask workers for parallel processing')
    
    parser.add_argument('--batch_size', type=int, default=200,
                       help='Number of time steps per processing batch. Reduce if OOM. (default: 200)')
    
    parser.add_argument('--compute_cloud_types', action='store_true', default=True,
                       help='Compute cloud type contributions')
    
    parser.add_argument('--skip_cloud_types', dest='compute_cloud_types',
                       action='store_false',
                       help='Skip cloud type computation')
    
    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_args()
    
    # Setup logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    print("="*80)
    print("SPATIAL STORM TYPE ATTRIBUTION FOR EXTREME PRECIPITATION")
    print("="*80)
    print(f"Catalog source: {args.catalog_source}")
    print(f"Percentiles: {', '.join(args.percentiles)}")
    print(f"Output directory: {args.output_dir}")
    print(f"Dask workers: {args.n_workers}")
    print(f"Batch size: {args.batch_size}")
    print(f"Compute cloud types: {args.compute_cloud_types}")
    print("="*80)
    
    # Load configuration
    print("\n📋 Loading configuration...")
    config = load_config(args.config_file, args.catalog_source)
    source_name = config.get('source_name')
    
    # Load storm masks
    root_dir = "/pscratch/sd/w/wcmca1/hackathon/"
    mask_dir = f"{root_dir}/cof_masks/{source_name}_cofmasks_hp8_v1.zarr"
    
    print(f"\n📂 Loading storm masks from: {mask_dir}")
    ds = xr.open_zarr(mask_dir, mask_and_scale=False, consolidated=True)
    ds = ds.pipe(egh.attach_coords)
    logger.info(f"Loaded mask dataset with {len(ds.time)} time steps")
    
    # Load extreme precipitation thresholds
    extreme_file = f"{root_dir}/extreme_precip/{source_name}_precip_percentiles_6h_hp8_v1.nc"
    print(f"📂 Loading extreme precipitation thresholds from: {extreme_file}")
    dsx = xr.open_dataset(extreme_file)
    
    # Load precipitation data
    catalog_file = "https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml"
    print(f"\n📂 Loading precipitation data...")
    pr = load_precipitation_data(catalog_file, args.catalog_source, config)
    logger.info(f"Loaded precipitation data with {len(pr.time)} time steps")
    
    # ========================================================================
    # CRITICAL: Align calendars and match time coordinates
    # ========================================================================
    logger.info("\n" + "="*60)
    logger.info("Checking and aligning time coordinates...")
    logger.info("="*60)
    
    # Standardize mask dataset calendar to match precipitation
    ds = standardize_calendar_type(ds, pr.time.values, logger)
    
    # Find common time range across datasets
    logger.info("\nFinding common time values between precipitation and mask datasets...")
    pr_times_set = set(pr.time.values)
    ds_times_set = set(ds.time.values)
    common_times = sorted(pr_times_set.intersection(ds_times_set))
    
    if not common_times:
        logger.error("❌ ERROR: No common time values between datasets!")
        logger.error(f"   Precipitation time range: {pr.time.values[0]} to {pr.time.values[-1]}")
        logger.error(f"   Mask time range: {ds.time.values[0]} to {ds.time.values[-1]}")
        return None
    
    logger.info(f"✅ Found {len(common_times)} common time points")
    logger.info(f"   Precipitation dataset: {len(pr.time)} time steps")
    logger.info(f"   Mask dataset: {len(ds.time)} time steps")
    logger.info(f"   Common times: {len(common_times)} time steps")
    
    # Select only the common times in both datasets
    pr = pr.sel(time=common_times)
    ds = ds.sel(time=common_times)
    logger.info("✅ Datasets aligned to common time values")
    
    # Apply time selection if specified
    if args.start_date and args.end_date:
        logger.info(f"\n⏱️  Subsetting time: {args.start_date} to {args.end_date}")
        pr = subset_time_range(pr.to_dataset(name='pr'), args.start_date, args.end_date, logger)['pr']
        ds = subset_time_range(ds, args.start_date, args.end_date, logger)
        date_suffix = f"_{args.start_date.replace('-', '')}_{args.end_date.replace('-', '')}"
        
        if len(pr.time) == 0 or len(ds.time) == 0:
            logger.error("❌ ERROR: No data found in specified time range!")
            return None
    else:
        # Get date range from data (for logging purposes)
        start_date_str = str(pr.time.values[0])[:10]
        end_date_str = str(pr.time.values[-1])[:10]
        args.start_date = start_date_str
        args.end_date = end_date_str
        date_suffix = ""  # No date suffix when using full dataset
    
    logger.info(f"✅ Final time alignment complete")
    logger.info(f"   Processing {len(pr.time)} time steps from {args.start_date} to {args.end_date}")
    logger.info("="*60 + "\n")
    
    # Verify dimensions match
    if len(pr.time) != len(ds.time):
        logger.error(f"❌ ERROR: Time dimension mismatch after alignment!")
        logger.error(f"   Precipitation: {len(pr.time)} time steps")
        logger.error(f"   Masks: {len(ds.time)} time steps")
        return None
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process each percentile
    for pname in args.percentiles:
        print(f"\n{'='*80}")
        print(f"PROCESSING PERCENTILE: {pname}")
        print(f"{'='*80}")
        
        # Get threshold
        threshold_var = f'pr_{pname.lower()}'
        if threshold_var not in dsx:
            print(f"⚠️  Warning: {threshold_var} not found in {extreme_file}")
            print(f"Available variables: {list(dsx.data_vars)}")
            continue
        
        pr_threshold = dsx[threshold_var].values
        print(f"Threshold value (median): {np.nanmedian(pr_threshold):.3f} mm/h")
        
        # Process with Dask
        start_time = time.time()
        ds_spatial = process_timeseries_dask(
            pr, pr_threshold, ds,
            percentile_name=pname,
            compute_cloud_types=args.compute_cloud_types,
            n_workers=args.n_workers,
            batch_size=args.batch_size
        )
        elapsed = time.time() - start_time
        print(f"\n⏱️  Processing time: {elapsed/60:.2f} minutes")
        
        # Save results
        output_file = f"{args.output_dir}/{source_name}_stormtype_spatial_{pname.lower()}{date_suffix}.nc"
        save_spatial_results(
            ds_spatial, output_file, source_name,
            args.start_date, args.end_date, pname
        )
    
    print("\n" + "="*80)
    print("✅ ALL PROCESSING COMPLETE!")
    print("="*80)


if __name__ == '__main__':
    main()
