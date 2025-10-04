#!/usr/bin/env python3
"""
Align NICAM Time Coordinates and Write Unified Dataset to Zarr
Only zoom 9 needs adjustment, as it has staggered time coordinates.

This script:
1. Loads NICAM data from intake catalog
2. Identifies variables with staggered time coordinates
3. Aligns all variables to a unified time coordinate
4. Writes the aligned dataset to a Zarr file using parallel dask

Variables included: pr, rlut, uas, vas

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import sys
import logging
import numpy as np
import pandas as pd
import xarray as xr
import intake
import easygems.healpix as egh
from pathlib import Path

# Add src directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
from zarr_tools import setup_dask_client, write_zarr


def setup_logging():
    """Configure logging for the script."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
        ]
    )
    return logging.getLogger(__name__)


def load_dataset(catalog_path, location, source, zoom, time):
    """
    Load dataset from intake catalog.
    
    Parameters:
    -----------
    catalog_path : str
        Path to intake catalog
    location : str
        Location key in catalog
    source : str
        Source identifier (e.g., 'nicam_gl11')
    zoom : int
        Zoom level
    time : str
        Time resolution (e.g., 'PT1H')
    
    Returns:
    --------
    xarray.Dataset
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Loading dataset: {source} from {location}")
    logger.info(f"  Zoom level: {zoom}, Time resolution: {time}")
    
    cat = intake.open_catalog(catalog_path)
    ds = cat[source](zoom=zoom, time=time).to_dask()
    ds = ds.pipe(egh.attach_coords)
    
    logger.info(f"  Dataset loaded: {len(ds.time)} timesteps, {len(ds.data_vars)} variables")
    logger.info(f"  Variables: {list(ds.data_vars.keys())}")
    
    return ds


def check_time_patterns(ds, logger=None):
    """
    Check which variables have data at :00 vs :30 timestamps.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Input dataset
    logger : logging.Logger, optional
        Logger for output
        
    Returns:
    --------
    tuple: (vars_at_30, vars_at_00)
    """
    if logger is None:
        logger = logging.getLogger(__name__)
        
    logger.info("Checking time patterns for all variables...")
    
    vars_at_30 = []
    vars_at_00 = []
    
    for var in ds.data_vars:
        # Check if variable has data at index 0 (:30) or index 1 (:00)
        count_at_0 = np.count_nonzero(~np.isnan(ds[var].isel(time=0).values))
        count_at_1 = np.count_nonzero(~np.isnan(ds[var].isel(time=1).values))
        
        if count_at_0 > 0:
            vars_at_30.append(var)
            logger.info(f"  {var:<10} → Valid at :30 timestamps (index 0, 2, 4, ...)")
        elif count_at_1 > 0:
            vars_at_00.append(var)
            logger.info(f"  {var:<10} → Valid at :00 timestamps (index 1, 3, 5, ...)")
        else:
            logger.warning(f"  {var:<10} → No valid data found in first timesteps!")
    
    logger.info(f"\nSummary:")
    logger.info(f"  Variables at :30: {vars_at_30}")
    logger.info(f"  Variables at :00: {vars_at_00}")
    
    return vars_at_30, vars_at_00


def create_unified_dataset(ds, vars_at_30min, vars_at_00min, logger=None):
    """
    Create a unified dataset with aligned time coordinates.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Original dataset with staggered times
    vars_at_30min : list
        Variables that have data at :30 timestamps (to be shifted)
    vars_at_00min : list
        Variables that have data at :00 timestamps (kept as-is)
    logger : logging.Logger, optional
        Logger for output
    
    Returns:
    --------
    xarray.Dataset with unified time coordinate
    
    Note:
    -----
    The first timestep of the original dataset (index 0) is typically :30,
    so we skip it to align with :00 variables that start at index 1.
    """
    if logger is None:
        logger = logging.getLogger(__name__)
        
    logger.info("Creating unified dataset with aligned time coordinates...")
    
    unified_vars = {}
    
    # Process variables at :30 (shift to :00)
    # Start at index 2 (not 0) to align with :00 variables at index 1
    logger.info(f"Processing {len(vars_at_30min)} variables at :30 (will shift by -30 min)...")
    for var in vars_at_30min:
        logger.info(f"  Processing {var}...")
        var_30 = ds[var].isel(time=slice(2, None, 2))
        new_time = var_30.time.values - pd.Timedelta(minutes=30)
        unified_vars[var] = var_30.assign_coords(time=new_time)
    
    # Process variables at :00 (keep as-is)
    logger.info(f"Processing {len(vars_at_00min)} variables at :00 (no shift needed)...")
    for var in vars_at_00min:
        logger.info(f"  Processing {var}...")
        unified_vars[var] = ds[var].isel(time=slice(1, None, 2))
    
    # Ensure all variables have the same length (trim to minimum)
    # This handles edge case where total timesteps might cause length mismatch
    min_len = min(len(unified_vars[var].time) for var in unified_vars)
    logger.info(f"Trimming all variables to minimum length: {min_len} timesteps")
    for var in unified_vars:
        unified_vars[var] = unified_vars[var].isel(time=slice(0, min_len))
    
    ds_unified = xr.Dataset(unified_vars)
    
    logger.info(f"Unified dataset created!")
    logger.info(f"  Original dataset: {len(ds.time)} timesteps")
    logger.info(f"  Unified dataset:  {len(ds_unified.time)} timesteps")
    logger.info(f"  Time range: {ds_unified.time.values[0]} to {ds_unified.time.values[-1]}")
    
    # Verify all variables have the same time coordinate
    all_times = [ds_unified[var].time.values for var in ds_unified.data_vars]
    all_match = all(np.array_equal(all_times[0], t) for t in all_times[1:])
    
    if all_match:
        logger.info("✓ All variables have matching time coordinates!")
    else:
        logger.error("✗ Time coordinates do not match across all variables!")
        raise ValueError("Time coordinate alignment failed!")
    
    return ds_unified


def main():
    """Main execution function."""
    logger = setup_logging()
    
    logger.info("="*80)
    logger.info("NICAM Time Coordinate Alignment Script")
    logger.info("="*80)
    
    # Configuration
    catalog_path = "/global/homes/f/feng045/program/hackathon/catalog/NERSC/main.yaml"
    location = "NERSC"
    source = "nicam_gl11"
    zoom = 9
    time = "PT1H"
    
    # Variables to include in unified dataset
    vars_to_include = ['pr', 'psl', 'rlut', 'rsut', 'uas', 'vas']
    
    # Output configuration
    output_zarr = f"/pscratch/sd/w/wcmca1/hackathon/NICAM/NICAM_2d1h_z{zoom}_shifted.zarr"
    
    # Dask configuration
    use_parallel = True
    n_workers = 8
    threads_per_worker = 4
    memory_per_worker = "30GB"
    
    # Chunk configuration for output
    chunksize_time = 24  # 24 hours per chunk
    chunksize_cell = None  # Let it use default
    
    logger.info("\nConfiguration:")
    logger.info(f"  Catalog: {catalog_path}")
    logger.info(f"  Source: {source}")
    logger.info(f"  Variables: {vars_to_include}")
    logger.info(f"  Output: {output_zarr}")
    logger.info(f"  Parallel: {use_parallel}")
    if use_parallel:
        logger.info(f"    Workers: {n_workers}")
        logger.info(f"    Threads per worker: {threads_per_worker}")
        logger.info(f"    Memory per worker: {memory_per_worker}")
    
    # Step 1: Load dataset
    logger.info("\n" + "="*80)
    logger.info("Step 1: Loading dataset")
    logger.info("="*80)
    ds = load_dataset(catalog_path, location, source, zoom, time)
    
    # Step 2: Check time patterns
    logger.info("\n" + "="*80)
    logger.info("Step 2: Checking time patterns")
    logger.info("="*80)
    vars_at_30, vars_at_00 = check_time_patterns(ds, logger)
    
    # Step 3: Filter to requested variables
    logger.info("\n" + "="*80)
    logger.info("Step 3: Filtering to requested variables")
    logger.info("="*80)
    vars_at_30_filtered = [v for v in vars_at_30 if v in vars_to_include]
    vars_at_00_filtered = [v for v in vars_at_00 if v in vars_to_include]
    
    logger.info(f"Variables to process:")
    logger.info(f"  At :30 (will shift): {vars_at_30_filtered}")
    logger.info(f"  At :00 (no shift):   {vars_at_00_filtered}")
    
    # Verify all requested variables are found
    all_filtered = vars_at_30_filtered + vars_at_00_filtered
    missing = set(vars_to_include) - set(all_filtered)
    if missing:
        logger.error(f"Missing variables: {missing}")
        raise ValueError(f"Requested variables not found in dataset: {missing}")
    
    # Step 4: Create unified dataset
    logger.info("\n" + "="*80)
    logger.info("Step 4: Creating unified dataset")
    logger.info("="*80)
    ds_unified = create_unified_dataset(ds, vars_at_30_filtered, vars_at_00_filtered, logger)
    
    # Step 5: Set up dask client
    logger.info("\n" + "="*80)
    logger.info("Step 5: Setting up Dask client")
    logger.info("="*80)
    client = setup_dask_client(
        parallel=use_parallel,
        n_workers=n_workers,
        threads_per_worker=threads_per_worker,
        memory_per_worker=memory_per_worker,
        logger=logger
    )
    
    if client:
        logger.info(f"Dask dashboard: {client.dashboard_link}")
    
    # Step 6: Write to Zarr
    logger.info("\n" + "="*80)
    logger.info("Step 6: Writing unified dataset to Zarr")
    logger.info("="*80)
    logger.info(f"Output path: {output_zarr}")
   
    try:
        write_zarr(
            ds_unified,
            output_zarr,
            client=client,
            logger=logger,
            chunksize_time=chunksize_time,
            chunksize_cell=chunksize_cell
        )
        logger.info("✓ Successfully wrote unified dataset to Zarr!")
        
    except Exception as e:
        logger.error(f"Error writing Zarr file: {e}")
        raise
    
    finally:
        # Clean up dask client
        if client:
            logger.info("Closing Dask client...")
            client.close()
    
    logger.info("\n" + "="*80)
    logger.info("COMPLETED SUCCESSFULLY")
    logger.info("="*80)


if __name__ == "__main__":
    main()
