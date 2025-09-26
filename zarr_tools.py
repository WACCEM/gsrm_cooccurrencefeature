#!/usr/bin/env python3
"""
Zarr Tools for Meteorological Data Processing

Common utilities for writing zarr files with optimized chunking and Dask support.
Shared functions for writing HEALPix-based meteorological datasets.

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import os
import shutil
import zarr
import numpy as np
import xarray as xr
import logging

try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

try:
    from dask.distributed import progress
    DASK_AVAILABLE = True
except ImportError:
    DASK_AVAILABLE = False


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
            # Temporarily suppress distributed.shuffle logs during progress display
            shuffle_logger = logging.getLogger('distributed.shuffle')
            original_level = shuffle_logger.level
            shuffle_logger.setLevel(logging.ERROR)  # Only show errors, not warnings
            
            # Get cluster state information before processing
            try:
                if PSUTIL_AVAILABLE:
                    memory_usage = client.run(lambda: psutil.Process().memory_info().rss / 1e9)
                    logger.info(f"Current memory usage across workers (GB): {memory_usage}")
            except Exception as e:
                logger.warning(f"Could not get memory usage: {e}")
                   
            try:
                # Compute with progress tracking
                future = client.compute(write_task)
                logger.info("Writing Zarr (this may take a while)...")
                if DASK_AVAILABLE:
                    progress(future)  # Shows a progress bar in notebooks or detailed progress in terminals

                result = future.result()
                logger.info("Zarr write completed successfully")
            except Exception as e:
                logger.error(f"Zarr write failed: {str(e)}")
                raise
            finally:
                # Restore original log level
                shuffle_logger.setLevel(original_level)
        except Exception as e:
            if not DASK_AVAILABLE:
                logger.warning("Dask distributed components not available, computing locally")
            else:
                logger.error(f"Error during distributed computation: {e}")
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
        logger.info(f"Removing existing zarr store: {out_zarr}")
        shutil.rmtree(out_zarr)
    
    logger.info(f"Writing zarr to: {out_zarr}")
    ds.to_zarr(out_zarr, mode='w', consolidated=True)
    logger.info(f"Zarr write completed: {out_zarr}")


def write_zarr_chunked(all_results, time_coords, mask_variables, output_path, 
                      template_coords, attrs, client=None, logger=None, chunk_size_time=24):
    """
    Write results to zarr in time chunks.
    
    Args:
        all_results: dict
            Dictionary with time keys and results
        time_coords: array-like
            Time coordinate values
        mask_variables: list
            List of mask variable names  
        output_path: str
            Path to output zarr file
        template_coords: dict
            Template coordinate dictionary
        attrs: dict
            Dataset attributes
        client: dask.distributed.Client, optional
            Dask client for distributed computation
        logger: logging.Logger, optional
            Logger for progress messages
        chunk_size_time: int, optional
            Number of time steps to process at once (default: 24)
    
    Returns:
        int: Number of successfully written time steps
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Remove existing zarr store if it exists
    if os.path.exists(output_path):
        logger.info(f"Removing existing zarr store: {output_path}")
        shutil.rmtree(output_path)
    
    # Ensure parent directory exists
    parent_dir = os.path.dirname(output_path)
    if parent_dir and not os.path.exists(parent_dir):
        os.makedirs(parent_dir, exist_ok=True)
        logger.info(f"Created parent directory: {parent_dir}")
    
    logger.info(f"Starting memory-efficient zarr write to: {output_path}")
    logger.info(f"Processing in chunks of {chunk_size_time} time steps")
    
    # Find a successful time step to get spatial dimensions and structure
    sample_time = None
    sample_results = None
    for time_val in time_coords:
        time_str = str(time_val)
        if time_str in all_results:
            sample_time = time_str
            sample_results = all_results[time_str]
            break
    
    if sample_results is None:
        logger.error("No successful results found in all_results")
        return 0
    
    logger.info(f"Using time step {sample_time} as template")
    
    # Get spatial dimensions from the first mask variable
    first_var = mask_variables[0]
    sample_data = sample_results[first_var]
    spatial_dims = sample_data.shape
    
    logger.info(f"Spatial dimensions: {spatial_dims}")
    logger.info(f"Total time steps to process: {len(time_coords)}")
    
    # Calculate optimal cell chunk size based on HEALPix zoom level
    chunksize_cell = None
    if 'crs' in template_coords:
        crs_coord = template_coords['crs']
        if hasattr(crs_coord, 'attrs') and 'healpix_nside' in crs_coord.attrs:
            zoom_level = zoom_level_from_nside(crs_coord.attrs['healpix_nside'])
            chunksize_cell = 12 * 4**zoom_level
            logger.info(f"HEALPix NSIDE={crs_coord.attrs['healpix_nside']}, zoom_level={zoom_level}, chunksize_cell={chunksize_cell}")
    
    if chunksize_cell is None:
        # Default cell chunking for non-HEALPix data
        chunksize_cell = min(10000, spatial_dims[0])
        logger.info(f"Using default cell chunking: {chunksize_cell}")
    
    # Log coordinate information for debugging
    if logger:
        logger.debug(f"Template coords keys: {list(template_coords.keys())}")
    
    # Handle CRS coordinate if present
    if 'crs' in template_coords:
        crs_coord = template_coords['crs']
        logger.debug(f"Found CRS coordinate: {crs_coord}")
    
    # Define the zarr store structure once
    zarr_store = None
    
    # Process time in chunks
    successful_times = 0
    total_chunks = (len(time_coords) + chunk_size_time - 1) // chunk_size_time
    
    logger.info(f"Will process {len(time_coords)} time steps in {total_chunks} chunks")
    
    # Main processing loop
    for chunk_idx in range(total_chunks):
        start_idx = chunk_idx * chunk_size_time
        end_idx = min(start_idx + chunk_size_time, len(time_coords))
        chunk_times = time_coords[start_idx:end_idx]
        
        logger.info(f"Processing chunk {chunk_idx + 1}/{total_chunks}: "
                   f"times {start_idx} to {end_idx-1}")
        
        # Collect data for this chunk
        chunk_data = {}
        chunk_time_coords = []
        
        try:
            for i, time_val in enumerate(chunk_times):
                time_str = str(time_val)
                if time_str in all_results:
                    timestep_results = all_results[time_str]
                    
                    # Initialize chunk_data structure on first successful time
                    if not chunk_data:
                        for var_name in mask_variables:
                            chunk_data[var_name] = []
                    
                    # Add data for each variable
                    for var_name in mask_variables:
                        var_data = timestep_results[var_name]
                        chunk_data[var_name].append(var_data)
                    
                    chunk_time_coords.append(time_val)
            
            if not chunk_data:
                logger.warning(f"No data found for chunk {chunk_idx + 1}, skipping")
                continue
            
            # Convert lists to numpy arrays with time dimension
            chunk_data_vars = {}
            for var_name in mask_variables:
                data_array = np.stack(chunk_data[var_name], axis=0)  # time is first dimension
                chunk_data_vars[var_name] = (['time', 'cell'], data_array)
                
            # Create other coordinates (excluding time which we'll set per chunk)
            other_coords = {k: v for k, v in template_coords.items()
                           if k not in ['time']}
            
            # Add time coordinate for this chunk
            chunk_coords = other_coords.copy()
            chunk_coords['time'] = chunk_time_coords
            
            # Create xarray Dataset for this chunk
            chunk_ds = xr.Dataset(
                data_vars=chunk_data_vars,
                coords=chunk_coords
            )
            
            # Apply chunking (optimized for HEALPix)
            chunk_ds = chunk_ds.chunk({'time': chunk_size_time, 'cell': chunksize_cell})
            
            # Set attributes
            if attrs:
                chunk_ds.attrs.update(attrs)
            
            # Write to zarr (append mode for subsequent chunks)
            if chunk_idx == 0:
                # First chunk - create zarr store
                logger.info("Creating new zarr store")
                chunk_ds.to_zarr(output_path, mode='w')
                zarr_store = output_path
                logger.info(f"Created zarr store with {len(chunk_time_coords)} time steps")
            else:
                # Subsequent chunks - append along time dimension  
                logger.info(f"Appending {len(chunk_time_coords)} time steps to zarr store")
                chunk_ds.to_zarr(output_path, mode='a', append_dim='time')
                logger.info(f"Appended chunk {chunk_idx + 1} successfully")
            
            successful_times += len(chunk_time_coords)
            
        except Exception as e:
            logger.error(f"Error processing chunk {chunk_idx + 1}: {e}")
            continue
            
        # Clear chunk data to free memory
        del chunk_data, chunk_data_vars, chunk_ds
        
        # Progress update
        if (chunk_idx + 1) % 5 == 0 or chunk_idx + 1 == total_chunks:
            logger.info(f"Completed {chunk_idx + 1}/{total_chunks} chunks")
    
    # Consolidate zarr metadata for better read performance
    try:
        zarr.consolidate_metadata(output_path)
        logger.info("Consolidated zarr metadata")
    except Exception as e:
        logger.warning(f"Could not consolidate zarr metadata: {e}")
    
    logger.info(f"Chunked zarr write completed: {successful_times} time steps written")
    return successful_times