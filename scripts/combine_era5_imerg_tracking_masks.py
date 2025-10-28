import xarray as xr
import pandas as pd
import numpy as np
import os
import glob
import time
import logging
import sys
from pathlib import Path
# Add src directory to path for zarr_tools import
sys.path.append(str(Path(__file__).parent.parent / 'src'))
from zarr_tools import setup_dask_client

#--------------------------------------------------------------------------------------
def setup_logging():
    """
    Set the logging message level

    Args:
        None.

    Returns:
        None.
    """
    logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO)

#--------------------------------------------------------------------------------------
def get_datasets(dir_mcs, dir_te, parallel=False, logger=None):
    """
    Load datasets from directories

    Args:
        dir_mcs: str
            Directory containing MCS Zarr files
        dir_te: str
            Directory containing TempestExtremes Zarr files
        parallel: bool
            Whether to use parallel processing
        logger: logging.Logger
            Logger for status messages
            
    Returns:
        tuple: (ds_mcs, ds_ar) datasets
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Read MCS file
    logger.info("Reading MCS file...")
    ds_mcs = xr.open_zarr(
        dir_mcs,
        consolidated=True,
        mask_and_scale=False,
    )
    logger.info(f"Finished reading MCS file")

    # Read TempestExtremes file
    logger.info("Reading TempestExtremes file...")
    ds_te = xr.open_zarr(
        dir_te,
        consolidated=True,
        mask_and_scale=False,
    )
    logger.info(f"Finished reading TempestExtremes file")

    return ds_mcs, ds_te

#--------------------------------------------------------------------------------------
def find_matching_times(target_times, source_times, tolerance):
    """Find times in source that match target within tolerance"""
    matching = []
    for t in target_times:
        # Find closest time in source
        time_diffs = np.abs(source_times - t)
        min_diff_idx = time_diffs.argmin()
        if time_diffs[min_diff_idx] <= tolerance:
            matching.append(source_times[min_diff_idx])
    return matching
    
#--------------------------------------------------------------------------------------
def combine_masks(ds_mcs, ds_te, client=None, out_zarr=None, logger=None):
    """
    Combine MCS and TempestExtremes tracking datasets.

    Args:
        ds_mcs: xarray.Dataset
            MCS tracking dataset
        ds_te: xarray.Dataset
            TempestExtremes tracking dataset
        client: dask.distributed.Client, optional
            Dask client for distributed computation
        out_zarr: str, optional
            Output Zarr store path
        logger: logging.Logger, optional
            Logger for status messages

    Returns:
        xarray.Dataset: Combined dataset with all masks
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    drop_var_list = ['']
    # rename_dict = {
    #     'AR_binary_tag': 'ar_mask',
    #     'TC_binary_tag': 'tc_mask',
    #     'ETC_binary_tag': 'etc_mask',
    #     'longitude': 'lon',
    #     'latitude': 'lat',
    # }
    
    # Find common time range across all datasets
    # Use pd.DatetimeIndex for robust time matching with tolerance    
    time_mcs = pd.DatetimeIndex(ds_mcs['time'].values)
    time_te = pd.DatetimeIndex(ds_te['time'].values)
    
    logger.info(f"MCS dataset has {len(time_mcs)} time steps")
    logger.info(f"TempestExtremes dataset has {len(time_te)} time steps")
    
    # Find overlapping time range
    start_time = max(time_mcs.min(), time_te.min())
    end_time = min(time_mcs.max(), time_te.max())
    
    # Create a common 6-hourly time coordinate from start to end
    # This ensures both datasets align to the same times
    common_times = pd.date_range(start=start_time, end=end_time, freq='6h')
    
    # Filter to only times that exist in both datasets (with small tolerance)
    tolerance = pd.Timedelta('10m')  # Allow 10 minutes difference for rounding
    
    # Find times present in both datasets
    mcs_matched = find_matching_times(common_times, time_mcs, tolerance)
    te_matched = find_matching_times(common_times, time_te, tolerance)
    
    # Get intersection of matched times (convert to timestamps for comparison)
    common_times = sorted(set(pd.DatetimeIndex(mcs_matched)).intersection(set(pd.DatetimeIndex(te_matched))))

    logger.info(f"Found {len(common_times)} common time steps")
    
    if not common_times:
        logger.warning("No common time values between all datasets!")
        return None
    else:
        # Select only the common times in all datasets (use nearest method for tolerance)
        ds_mcs = ds_mcs.sel(time=common_times, method='nearest', tolerance=tolerance)
        ds_te = ds_te.sel(time=common_times, method='nearest', tolerance=tolerance)

    # Fix for lat/lon coordinates issue: ensure consistent treatment
    datasets = [ds_mcs, ds_te]

    # Merge the datasets
    ds = xr.merge(datasets, combine_attrs='drop_conflicts', compat='override')
    logger.info(f"Successfully merged datasets with {len(common_times)} common time points")

    # Rename variables, drop unwanted ones in the DataSet
    ds = ds.drop_vars(drop_var_list, errors='ignore')
    # ds = ds.rename(rename_dict).drop_vars(drop_var_list, errors='ignore')

    # Remove _FillValue from attributes and set it in encoding instead
    # This prevents conflicts when xarray tries to encode the variables
    for var in ds.data_vars:
        # Remove _FillValue from attributes if it exists
        if '_FillValue' in ds[var].attrs:
            fill_value = ds[var].attrs.pop('_FillValue')
        else:
            fill_value = 0  # Default for mask variables
        
        # Set fill value in encoding based on dtype
        if np.issubdtype(ds[var].dtype, np.integer):
            ds[var].encoding['_FillValue'] = int(fill_value)
        elif np.issubdtype(ds[var].dtype, np.floating):
            # Use 0.0 for float mask variables instead of NaN
            ds[var].encoding['_FillValue'] = 0.0 if np.isnan(fill_value) else float(fill_value)
        
        # Clear other encoding that might cause conflicts
        for key in ['chunks', 'preferred_chunks']:
            ds[var].encoding.pop(key, None)
    
    # Also handle coordinates
    for coord in ds.coords:
        if coord in ds.variables:
            # Remove _FillValue from coordinate attributes
            if '_FillValue' in ds[coord].attrs:
                ds[coord].attrs.pop('_FillValue')
            # Clear encoding that might cause conflicts
            for key in ['chunks', 'preferred_chunks', '_FillValue']:
                ds[coord].encoding.pop(key, None)

    # # Rechunk to ensure consistent chunking across all variables
    # # Fix the inconsistent cell chunking that results from the merge
    # # Get the HEALPix zoom level to calculate proper cell chunk size
    # zoom_level = zoom_level_from_nside(ds.crs.attrs['healpix_nside'])
    # chunksize_cell = 12 * 4**zoom_level
    
    # logger.info(f"Before rechunk: {dict(ds.chunks)}")
    # ds = ds.chunk({'time': 28, 'cell': chunksize_cell})
    # logger.info(f"After rechunk: {dict(ds.chunks)}")

    # Add global attributes
    ds.attrs['processing_date'] = str(np.datetime64('today'))
    ds.attrs['processing_script'] = os.path.basename(__file__)

    # Write to Zarr
    if out_zarr:
        write_zarr(ds, out_zarr, client=client, logger=logger)

    return ds


#-------------------------------------------------------------------
def write_zarr(ds, out_zarr, client=None, logger=None):
    """
    Write dataset to Zarr with optimized chunking for HEALPix grid.
    
    Args:
        ds: xarray.Dataset
            Dataset to write
        out_zarr: str
            Output Zarr store path
        client: dask.distributed.Client, optional
            Dask client for distributed computation
        logger: logging.Logger, optional
            Logger for status messages
            
    Returns:
        None
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Optimize cell chunking for HEALPix grid
    zoom_level = zoom_level_from_nside(ds.crs.attrs['healpix_nside'])
    chunksize_time = 28
    chunksize_cell = 12 * 4**zoom_level
    
    # Make time chunks more even if needed
    if isinstance(chunksize_time, (int, float)) and chunksize_time != 'auto':
        total_times = ds.sizes['time']
        chunks = total_times // chunksize_time
        if chunks * chunksize_time < total_times:
            # We have a remainder - try to make chunks more even
            if total_times % chunks == 0:
                chunksize_time = total_times // chunks
            elif total_times % (chunks + 1) == 0:
                chunksize_time = total_times // (chunks + 1)
    
    # Set proper chunking for HEALPix output
    chunked_hp = ds.chunk({
        "time": chunksize_time, 
        "cell": chunksize_cell, 
    })
    # Report dataset size and chunking info
    logger.info(f"Output dataset dimensions: {dict(chunked_hp.sizes)}")
    logger.info(f"Output chunking scheme: time={chunksize_time}, cell={chunksize_cell}")

    # ---------- WRITE HEALPIX ZARR OUTPUT ----------
    logger.info(f"Starting Zarr write to: {out_zarr}")
    
    # Create a delayed task for Zarr writing
    write_task = chunked_hp.to_zarr(
        out_zarr,
        mode="w",
        consolidated=True,  # Enable for better performance when reading
        compute=False      # Create a delayed task
    )
    
    # Compute the task, with progress reporting
    if client:
        from dask.distributed import progress
        import psutil

        # Temporarily suppress distributed.shuffle logs during progress display
        shuffle_logger = logging.getLogger('distributed.shuffle')
        original_level = shuffle_logger.level
        shuffle_logger.setLevel(logging.ERROR)  # Only show errors, not warnings

        # Get cluster state information before processing
        memory_usage = client.run(lambda: psutil.Process().memory_info().rss / 1e9)
        logger.info(f"Current memory usage across workers (GB): {memory_usage}")
               
        try:
            # Compute with progress tracking
            future = client.compute(write_task)
            logger.info("Writing Zarr (this may take a while)...")
            progress(future)  # Shows a progress bar in notebooks or detailed progress in terminals

            result = future.result()
            logger.info("Zarr write completed successfully")
        except Exception as e:
            logger.error(f"Zarr write failed: {str(e)}")
            raise
        finally:
            # Restore original log level
            shuffle_logger.setLevel(original_level)
    else:
        # Compute locally if no client
        write_task.compute()

    logger.info(f"✅ Zarr file complete: {out_zarr}")

#-------------------------------------------------------------------
def zoom_level_from_nside(nside):
    """
    Calculate the zoom level from the NSIDE value.

    Args:
        nside (int): NSIDE value, must be a power of 2.
    
    Returns:
        int: Zoom level corresponding to the NSIDE value.
    """
    zoom = int(np.log2(nside))
    if 2**zoom != nside:
        raise ValueError("NSIDE must be a power of 2.")
    return zoom


def main():
    """Main function to run the remap masks"""
    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)

    start_time = time.time()
    logger.info("Starting remap masks ...")

    # Define parameters
    source_name = "IMERGv7"
    zoom = 8
    version = "v1"
    parallel = True
    n_workers = 8
    threads_per_worker = 4

    in_dir = "/pscratch/sd/w/wcmca1/hackathon/all_masks/"
    dir_te = f"{in_dir}ERA5_AR_TC_ETC_hp{zoom}_{version}.zarr"
    dir_mcs = f"/pscratch/sd/w/wcmca1/hackathon/mcs_masks/{source_name}_mcs_masks_hp{zoom}.zarr"

    out_dir = "/pscratch/sd/w/wcmca1/hackathon/all_masks/"
    out_basename = f"{source_name}_allmasks_hp{zoom}_{version}.zarr"
    out_zarr = f"{out_dir}{out_basename}"
    os.makedirs(out_dir, exist_ok=True)

    # Setup Dask client
    client = setup_dask_client(parallel=parallel, n_workers=n_workers, threads_per_worker=threads_per_worker, logger=logger)

    try:
        # Load datasets
        ds_mcs, ds_te = get_datasets(dir_mcs, dir_te, parallel, logger)

        # Combine datasets
        ds = combine_masks(ds_mcs, ds_te, client=client, out_zarr=out_zarr, logger=logger)
        
        # Cleanup
        ds_mcs.close()
        ds_te.close()
        ds.close()
        
    finally:
        # Always cleanup client
        if client and parallel:
            # Suppress Dask shutdown messages by temporarily raising log level
            logging.getLogger('distributed').setLevel(logging.CRITICAL)
            logging.getLogger('distributed.worker').setLevel(logging.CRITICAL)
            logging.getLogger('distributed.nanny').setLevel(logging.CRITICAL)
            
            client.close()
    
    # Log completion time
    end_time = time.time()
    elapsed_time = end_time - start_time
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    print(f"\n{'='*80}")
    print(f"✅ PROCESSING COMPLETE!")
    print(f"{'='*80}")
    print(f"Output: '{out_zarr}'")
    print(f"Total time: {int(seconds):02} seconds ({int(minutes):02} minutes)")

if __name__ == "__main__":
    main()