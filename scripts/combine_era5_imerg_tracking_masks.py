import xarray as xr
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
import intake
import requests
import easygems.healpix as egh
from functools import partial

#-------------------------------------------------------------------
def setup_logging():
    """
    Set the logging message level

    Args:
        None.

    Returns:
        None.
    """
    logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO)

#-------------------------------------------------------------------
# Dask client setup moved to zarr_tools.py - import it instead

#-------------------------------------------------------------------
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

#-------------------------------------------------------------------
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

    drop_var_list = ['ccs_mask']
    # rename_dict = {
    #     'AR_binary_tag': 'ar_mask',
    #     'TC_binary_tag': 'tc_mask',
    #     'ETC_binary_tag': 'etc_mask',
    #     'longitude': 'lon',
    #     'latitude': 'lat',
    # }
    
    # Find common time range across all datasets
    common_times = sorted(set(ds_mcs['time'].values)
                         .intersection(set(ds_te['time'].values)))
    if not common_times:
        logger.warning("No common time values between all datasets!")
        return None
    else:
        # Select only the common times in all datasets
        ds_mcs = ds_mcs.sel(time=common_times)
        ds_te = ds_te.sel(time=common_times)

    # Fix for lat/lon coordinates issue: ensure consistent treatment
    datasets = [ds_mcs, ds_te]

    # Merge the datasets
    ds = xr.merge(datasets, combine_attrs='drop_conflicts', compat='override')
    logger.info(f"Successfully merged datasets with {len(common_times)} common time points")

    # Rename variables, drop unwanted ones in the DataSet
    ds = ds.drop_vars(drop_var_list, errors='ignore')
    # ds = ds.rename(rename_dict).drop_vars(drop_var_list, errors='ignore')

    # TODO: Modify global attributes if needed
    # ds.attrs['history'] = f"Created on {time.ctime()} by combining tracking data"

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
    chunksize_time = 24
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

    logger.info(f"Zarr file complete: {out_zarr}")

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
    threads_per_worker = 16
    memory_per_worker = "60GB"
    # chunksize_cell = 12 * 4**zoom
    # chunksize_time = 24

    # catalog_dict = {
    #     "catalog_file": "https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml",
    #     "catalog_location": "NERSC",
    #     "catalog_source": "IR_IMERG",
    #     "catalog_params": {"zoom": zoom},
    # }

    in_dir = "/pscratch/sd/w/wcmca1/hackathon/allmasks/"
    dir_te = f"{in_dir}ERA5_AR_TC_ETC_hp{zoom}_{version}.zarr"
    dir_mcs = f"/pscratch/sd/w/wcmca1/hackathon/mcs/{source_name}/mcstracking/{source_name}_hrly_mcsmask_hp{zoom}_v1.zarr"

    out_dir = "/pscratch/sd/w/wcmca1/hackathon/allmasks/"
    out_basename = f"{source_name}_allmasks_hp{zoom}_{version}.zarr"
    out_zarr = f"{out_dir}{out_basename}"
    os.makedirs(out_dir, exist_ok=True)

    # Setup Dask client
    client = setup_dask_client(parallel, n_workers, threads_per_worker, memory_per_worker, logger)

    try:
        # Load datasets
        ds_mcs, ds_te = get_datasets(dir_mcs, dir_te, parallel, logger)

        # Combine datasets
        ds = combine_masks(ds_mcs, ds_te, client=client, out_zarr=out_zarr, logger=logger)
        # import pdb; pdb.set_trace()
        
        # Cleanup
        ds_mcs.close()
        ds_te.close()
        ds.close()
        
    finally:
        # Always cleanup client
        if client and parallel:
            logger.info("Shutting down Dask client")
            client.close()

    # Log completion time
    end_time = time.time()
    elapsed_time = end_time - start_time
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    logger.info(f"Conversion completed in {int(hours):02}:{int(minutes):02}:{int(seconds):02} (hh:mm:ss).")

if __name__ == "__main__":
    main()