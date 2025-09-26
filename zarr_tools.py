#!/usr/bin/env python3
"""
Zarr Tools for Meteorological Data Processing

Common utilities for writing zarr files with optimized chunking and Dask support.
Shared functions for writing HEALPix-based meteorological datasets.

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import os
import numpy as np
import xarray as xr
import logging


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


def write_zarr(ds, out_zarr, client=None, logger=None, chunksize_time=24, chunksize_cell=None):
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
        chunksize_time: int, optional
            Time chunk size (default: 24)
        chunksize_cell: int, optional
            Cell chunk size. If None, calculated from HEALPix nside
            
    Returns:
        None
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Remove existing zarr store if it exists
    if os.path.exists(out_zarr):
        import shutil
        logger.info(f"Removing existing zarr store: {out_zarr}")
        shutil.rmtree(out_zarr)

    # Optimize cell chunking for HEALPix grid if not provided
    if chunksize_cell is None:
        if hasattr(ds, 'crs') and 'healpix_nside' in ds.crs.attrs:
            zoom_level = zoom_level_from_nside(ds.crs.attrs['healpix_nside'])
            chunksize_cell = 12 * 4**zoom_level
        else:
            # Default cell chunking for non-HEALPix data
            chunksize_cell = min(10000, ds.sizes.get('cell', ds.sizes.get('ncol', 1000)))
    
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
    
    # Determine spatial dimension name
    spatial_dim = None
    for dim_name in ['cell', 'ncol']:
        if dim_name in ds.dims:
            spatial_dim = dim_name
            break
    
    if spatial_dim is None:
        raise ValueError("Could not find spatial dimension ('cell' or 'ncol') in dataset")
    
    # Set proper chunking
    chunk_dict = {
        "time": chunksize_time, 
        spatial_dim: chunksize_cell, 
    }
    chunked_ds = ds.chunk(chunk_dict)
    
    # Report dataset size and chunking info
    logger.info(f"Output dataset dimensions: {dict(chunked_ds.sizes)}")
    logger.info(f"Output chunking scheme: time={chunksize_time}, {spatial_dim}={chunksize_cell}")

    # ---------- WRITE ZARR OUTPUT ----------
    logger.info(f"Starting Zarr write to: {out_zarr}")
    
    # Create a delayed task for Zarr writing
    write_task = chunked_ds.to_zarr(
        out_zarr,
        mode="w",
        consolidated=True,  # Enable for better performance when reading
        compute=False      # Create a delayed task
    )
    
    # Compute the task, with progress reporting
    if client:
        try:
            from dask.distributed import progress
            import psutil

            # Temporarily suppress distributed.shuffle logs during progress display
            shuffle_logger = logging.getLogger('distributed.shuffle')
            original_level = shuffle_logger.level
            shuffle_logger.setLevel(logging.ERROR)  # Only show errors, not warnings

            # Get cluster state information before processing
            try:
                memory_usage = client.run(lambda: psutil.Process().memory_info().rss / 1e9)
                logger.info(f"Current memory usage across workers (GB): {memory_usage}")
            except Exception as e:
                logger.warning(f"Could not get memory usage: {e}")
                   
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
        except ImportError:
            logger.warning("Dask distributed components not available, computing locally")
            write_task.compute()
    else:
        # Compute locally if no client
        logger.info("Computing zarr write locally...")
        write_task.compute()

    logger.info(f"Zarr file complete: {out_zarr}")


def write_zarr_simple(ds, out_zarr, logger=None):
    """
    Simple zarr write function without Dask optimization.
    
    Args:
        ds: xarray.Dataset
            Dataset to write
        out_zarr: str
            Output Zarr store path
        logger: logging.Logger, optional
            Logger for status messages
            
    Returns:
        None
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Remove existing zarr store if it exists
    if os.path.exists(out_zarr):
        import shutil
        logger.info(f"Removing existing zarr store: {out_zarr}")
        shutil.rmtree(out_zarr)
    
    logger.info(f"Writing zarr to: {out_zarr}")
    ds.to_zarr(out_zarr, mode='w', consolidated=True)
    logger.info(f"Zarr write completed: {out_zarr}")