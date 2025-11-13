#!/usr/bin/env python3
"""
Zarr Tools for Meteorological Data Processing

Common utilities for writing zarr files with optimized chunking and Dask support.
Shared functions for writing HEALPix-based meteorological datasets.

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import os
import shutil
import sys
import traceback
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
    from dask.distributed import progress, as_completed
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
    
    # Clear any existing encoding that might conflict with new chunking
    # This is critical when rechunking data that was previously written to zarr
    for var in ds.data_vars:
        if 'chunks' in ds[var].encoding:
            del ds[var].encoding['chunks']
        if 'preferred_chunks' in ds[var].encoding:
            del ds[var].encoding['preferred_chunks']
    
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


def initialize_zarr_store(output_path, time_coords, mask_variables, template_coords, attrs, 
                         chunk_size_time=24, var_attrs=None):
    """
    Initialize zarr store with proper structure for streaming writes using xarray.
    
    This approach creates a proper xarray Dataset first, then saves to zarr,
    which ensures coordinates are handled correctly.
    
    Parameters:
    -----------
    output_path : str
        Path where zarr store will be created
    time_coords : numpy.ndarray
        All time coordinates for the dataset
    mask_variables : list
        List of mask variable names to create
    template_coords : dict
        Coordinate information from template dataset
    attrs : dict
        Global attributes for the dataset
    chunk_size_time : int
        Time dimension chunk size
    var_attrs : dict, optional
        Dictionary mapping variable names to their attributes
        Example: {'mcs_mask': {'long_name': 'MCS swath mask', 'units': '1'}}
    """
    
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    
    # Get spatial dimensions from template
    n_cells = len(template_coords['cell'])
    
    # Calculate optimal cell chunking for HEALPix
    chunksize_cell = None
    if 'crs' in template_coords and hasattr(template_coords['crs'], 'attrs') and 'healpix_nside' in template_coords['crs'].attrs:
        zoom_level = zoom_level_from_nside(template_coords['crs'].attrs['healpix_nside'])
        chunksize_cell = 12 * 4**zoom_level
    else:
        # Default cell chunking for non-HEALPix data
        chunksize_cell = min(10000, n_cells)
    
    # Create data variables dictionary
    data_vars = {}
    for var_name in mask_variables:
        # Start with default attributes
        default_attrs = {'grid_mapping': 'crs', '_FillValue': 0.0}
        
        # Add user-specified attributes if provided
        if var_attrs and var_name in var_attrs:
            default_attrs.update(var_attrs[var_name])
        
        # Initialize with zeros
        data_vars[var_name] = (['time', 'cell'], 
                              np.zeros((len(time_coords), n_cells), dtype=np.float32),
                              default_attrs)
    
    # Create coordinates dictionary (only the actual coordinates)
    coords = {
        'time': (('time',), time_coords, dict(template_coords['time'].attrs)),
        'cell': (('cell',), template_coords['cell'].values.astype(int), dict(template_coords['cell'].attrs))
    }
    
    # Add crs coordinate if it exists
    if 'crs' in template_coords:
        coords['crs'] = (template_coords['crs'].dims, 
                        template_coords['crs'].values, 
                        dict(template_coords['crs'].attrs))
    
    # Create xarray Dataset
    ds = xr.Dataset(
        data_vars=data_vars,
        coords=coords,
        attrs=attrs
    )
    
    # Define chunking for zarr and chunk the dataset
    chunks = {
        'time': chunk_size_time,
        'cell': chunksize_cell
    }
    ds = ds.chunk(chunks)
    
    # Write to zarr
    ds.to_zarr(output_path, mode='w')
    
    print(f"Initialized zarr store: {len(time_coords)} time steps, {n_cells} cells")
    
    # Close the dataset to free memory
    ds.close()


def append_chunk_to_zarr(chunk_results, chunk_times, chunk_idx, mask_variables, output_path, logger=None):
    """
    Append a chunk of processed results to the existing zarr store.
    
    Parameters:
    -----------
    chunk_results : dict
        Dictionary mapping time strings to result dictionaries
    chunk_times : numpy.ndarray
        Time coordinates for this chunk
    chunk_idx : int
        Index of this chunk (for determining time slice)
    mask_variables : list
        List of mask variable names to write
    output_path : str
        Path to zarr store
    logger : logging.Logger, optional
        Logger for status messages
    """
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Open existing zarr store
    zarr_store = zarr.open(output_path, mode='r+')
    
    # Calculate time indices for this chunk
    time_start = chunk_idx * len(chunk_times)  # Assumes uniform chunk sizes
    time_end = time_start + len(chunk_times)
    
    # Write each mask variable
    for var_name in mask_variables:
        var_data = []
        
        # Collect data for this variable across all time steps in chunk
        for time_val in chunk_times:
            time_str = str(time_val)
            if time_str in chunk_results and var_name in chunk_results[time_str]:
                # Get the mask data
                mask_data = chunk_results[time_str][var_name]
                if hasattr(mask_data, 'values'):
                    mask_data = mask_data.values
                var_data.append(mask_data)
            else:
                # Fill with zeros if data missing
                n_cells = zarr_store[var_name].shape[1]
                var_data.append(np.zeros(n_cells, dtype=np.float32))
        
        # Convert to numpy array and write to zarr
        if var_data:
            var_array = np.stack(var_data, axis=0)
            zarr_store[var_name][time_start:time_end, :] = var_array
    
    logger.info(f"Appended chunk {chunk_idx + 1} data to zarr store")


def stream_process_to_zarr(ds, time_coords, mask_variables, output_path, template_coords, attrs,
                          client=None, logger=None, parallel=True, chunk_size_time=24, input_zarr_path=None):
    """
    Stream processing and writing to zarr without accumulating all results in memory.
    
    This function processes time steps in chunks, writes each chunk to zarr immediately,
    and frees memory before processing the next chunk. This eliminates the memory
    accumulation issue where all_results grows to 157GB.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Full input dataset
    time_coords : numpy.ndarray
        Time coordinates to process
    mask_variables : list
        List of mask variable names to create
    output_path : str
        Path for zarr output
    template_coords : dict
        Coordinate information from template dataset
    attrs : dict
        Global attributes for the output dataset
    client : dask.distributed.Client, optional
        Dask client for parallel processing
    logger : logging.Logger, optional
        Logger for status messages
    parallel : bool
        Whether to use parallel processing
    chunk_size_time : int
        Number of time steps to process in each chunk
        
    Returns:
    --------
    tuple : (int, list) - Number of successfully processed time steps and list of ETC overlap records
    """
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Process and write time steps in chunks
    total_processed = 0
    total_chunks = (len(time_coords) + chunk_size_time - 1) // chunk_size_time
    all_etc_records = []  # Collect ETC overlap records from all timesteps
    
    logger.info(f"Processing {len(time_coords)} time steps in {total_chunks} chunks of {chunk_size_time}")
    
    for chunk_idx in range(total_chunks):
        start_idx = chunk_idx * chunk_size_time
        end_idx = min((chunk_idx + 1) * chunk_size_time, len(time_coords))
        chunk_times = time_coords[start_idx:end_idx]
        
        logger.info(f"Processing chunk {chunk_idx + 1}/{total_chunks}: time steps {start_idx}-{end_idx-1}")
        
        # Process this chunk of time steps
        chunk_results = {}
        
        if parallel and client is not None:
            # Add current directory to path to import the main script functions
            current_dir = os.path.dirname(os.path.abspath(__file__))
            if current_dir not in sys.path:
                sys.path.insert(0, current_dir)            
        
        if parallel and client is not None:
            # Parallel processing for this chunk using zarr file path approach
            futures = []
            serialization_failed = False
            
            logger.info(f"Submitting {len(chunk_times)} tasks to {len(client.nthreads())} workers...")
            
            # Use the provided input zarr path
            if input_zarr_path is None:
                logger.error("No input zarr path provided. Falling back to sequential processing.")
                serialization_failed = True
            else:
                logger.info(f"Using input zarr file path: {input_zarr_path}")
                
                try:
                    for time_val in chunk_times:
                        # Submit with just time_val and zarr_path - minimal serialization!
                        future = client.submit(process_timestep_wrapper_zarr, time_val, input_zarr_path, verbose=False)
                        futures.append(future)
                        
                    logger.info(f"Successfully submitted {len(futures)} tasks with minimal serialization")
                        
                except Exception as e:
                    logger.error(f"Error submitting zarr-based tasks: {e}")
                    serialization_failed = True
                    futures = []
            
            if not serialization_failed and futures:
                # Collect results for this chunk
                try:
                    logger.info(f"Collecting results from {len(futures)} parallel tasks...")
                    
                    for future in as_completed(futures):
                        try:
                            time_str, timestep_results = future.result()
                            if timestep_results is not None:
                                chunk_results[time_str] = timestep_results
                                # Extract ETC overlap records if present
                                if 'etc_overlap_records' in timestep_results:
                                    all_etc_records.extend(timestep_results['etc_overlap_records'])
                            else:
                                logger.warning(f"Failed to process time step: {time_str}")
                        except Exception as e:
                            logger.error(f"Error processing time step: {e}")
                            # Continue processing other time steps rather than failing the entire chunk
                            
                except ImportError:
                    logger.error("Dask as_completed not available. This should not happen if we got this far.")
                    serialization_failed = True
                    
            # Fall back to sequential processing only if input zarr path missing or dask unavailable        
            if serialization_failed or not parallel or client is None:
                # Sequential processing for this chunk
                logger.info(f"Processing chunk {chunk_idx + 1} sequentially...")
                try:
                    from make_cooccurrence_masks import process_single_timestep_overlaps
                except ImportError:
                    logger.error("Could not import processing function. Aborting.")
                    return 0
                    
                for time_val in chunk_times:
                    try:
                        _ds = ds.sel(time=time_val)
                        timestep_results = process_single_timestep_overlaps(_ds, verbose=False)
                        chunk_results[str(time_val)] = timestep_results
                        # Extract ETC overlap records if present
                        if 'etc_overlap_records' in timestep_results:
                            all_etc_records.extend(timestep_results['etc_overlap_records'])
                    except Exception as e:
                        logger.error(f"Error processing time step {time_val}: {e}")
        
        # Write this chunk to zarr immediately
        try:
            logger.info(f"Writing chunk {chunk_idx + 1} to zarr...")
            append_chunk_to_zarr(
                chunk_results=chunk_results,
                chunk_times=chunk_times,
                chunk_idx=chunk_idx,
                mask_variables=mask_variables,
                output_path=output_path,
                logger=logger
            )
            
            # Update progress
            processed_this_chunk = len(chunk_results)
            total_processed += processed_this_chunk
            logger.info(f"Chunk {chunk_idx + 1} complete: {processed_this_chunk}/{len(chunk_times)} time steps written")
            
            # Free memory by explicitly deleting chunk results
            del chunk_results
            
        except Exception as e:
            logger.error(f"Error writing chunk {chunk_idx + 1} to zarr: {e}")
            continue
        
        # Log memory usage periodically
        if (chunk_idx + 1) % 5 == 0:
            if PSUTIL_AVAILABLE:
                try:
                    memory_usage = psutil.virtual_memory().percent
                    logger.info(f"Memory usage after chunk {chunk_idx + 1}: {memory_usage:.1f}%")
                except:
                    pass
    
    # Final consolidation of metadata
    try:
        zarr.consolidate_metadata(output_path)
        logger.info("Consolidated zarr metadata")
    except Exception as e:
        logger.warning(f"Could not consolidate zarr metadata: {e}")
    
    logger.info(f"Streaming processing complete: {total_processed}/{len(time_coords)} time steps successful")
    logger.info(f"Collected {len(all_etc_records)} ETC overlap records")
    return total_processed, all_etc_records


def process_timestep_wrapper_zarr(time_val, zarr_path, verbose=False):
    """
    Wrapper function for processing a single time step by reading from zarr file.
    
    This approach avoids serialization issues by having each worker read the zarr file
    directly instead of serializing large xarray Datasets.
    
    Args:
        time_val: Time coordinate value to process
        zarr_path: Path to the zarr file containing the data
        verbose: Whether to print verbose output
        
    Returns:
        tuple: (time_str, results_dict) or (time_str, None) if error
    """
    try:
        # Import here to avoid circular imports
        from make_cooccurrence_masks import process_single_timestep_overlaps
        
        # Each worker opens the zarr file independently
        ds = xr.open_dataset(zarr_path, engine='zarr')
        
        # Select the specific time step
        _ds = ds.sel(time=time_val)
        
        # Load the data for this time step into memory
        _ds = _ds.load()
        
        # Close the full dataset to free memory
        ds.close()
        
        # Process this time step
        timestep_results = process_single_timestep_overlaps(_ds, verbose=verbose)
        
        return str(time_val), timestep_results
        
    except Exception as e:
        print(f"Error processing time step {time_val}: {e}")
        traceback.print_exc()
        return str(time_val), None

def setup_dask_client(parallel, n_workers, threads_per_worker, memory_per_worker=None, 
                      dashboard_address=':8787', scheduler_port=0, processes=True,
                      memory_target_fraction=0.8, memory_spill_fraction=0.85, 
                      memory_pause_fraction=0.9, memory_terminate_fraction=0.95,
                      tcp_timeout='300s', client_heartbeat='10s', logger=None):
    """
    Set up a Dask client optimized for HPC hardware with comprehensive configuration options
    
    Args:
        parallel: bool
            Whether to use parallel processing
        n_workers: int
            Number of workers for the Dask cluster
        threads_per_worker: int
            Number of threads per worker
        memory_per_worker: str, optional
            Memory limit per worker (e.g., "60GB"). If None, auto-calculated from system memory
        dashboard_address: str, optional
            Dashboard address (default: ':8787')
        scheduler_port: int, optional
            Scheduler port (default: 0 for auto-selection)
        processes: bool, optional
            Use processes for better memory isolation (default: True)
        memory_target_fraction: float, optional
            Target memory fraction before spilling to disk (default: 0.8)
        memory_spill_fraction: float, optional
            Memory fraction at which to start spilling to disk (default: 0.85)
        memory_pause_fraction: float, optional
            Memory fraction at which to pause worker (default: 0.9)
        memory_terminate_fraction: float, optional
            Memory fraction at which to terminate worker (default: 0.95)
        tcp_timeout: str, optional
            TCP timeout for distributed communication (default: '300s')
        client_heartbeat: str, optional
            Client heartbeat interval (default: '10s')
        logger: logging.Logger, optional
            Logger for status messages
            
    Returns:
        dask.distributed.Client or None: Dask client if parallel is True, None otherwise
    """
    if logger is None:
        logger = logging.getLogger(__name__)
        
    if not parallel:
        logger.info("Running in sequential mode (parallel=False)")
        return None
    
    try:
        from dask.distributed import Client, LocalCluster
        import dask
        import psutil
        
        # Auto-calculate memory limit per worker if not provided
        if memory_per_worker is None:
            total_memory_gb = psutil.virtual_memory().total / (1024**3)  # Convert to GB
            # Leave 20% for system overhead, divide remaining by number of workers
            usable_memory_gb = total_memory_gb * 0.8
            memory_per_worker_gb = usable_memory_gb / n_workers
            memory_per_worker = f"{memory_per_worker_gb:.1f}GB"
            logger.info(f"Auto-calculated memory per worker: {memory_per_worker} "
                       f"(from {total_memory_gb:.1f}GB total system memory)")
        
        logger.info(f"Setting up Dask cluster optimized for HPC hardware")
        logger.info(f"Workers: {n_workers}, Threads per worker: {threads_per_worker}")
        logger.info(f"Memory per worker: {memory_per_worker}")
        logger.info(f"Memory management: target={memory_target_fraction}, "
                   f"spill={memory_spill_fraction}, pause={memory_pause_fraction}")
        
        # Set up environment variables for NUMA-aware computation (commented out by default)
        # import os
        # os.environ['OMP_NUM_THREADS'] = str(threads_per_worker)
        # os.environ['MKL_NUM_THREADS'] = str(threads_per_worker)
        # os.environ['OPENBLAS_NUM_THREADS'] = str(threads_per_worker)
        # os.environ['NUMBA_NUM_THREADS'] = str(threads_per_worker)
        
        # Configure Dask settings before creating cluster (new approach)
        dask.config.set({
            'distributed.worker.memory.target': memory_target_fraction,
            'distributed.worker.memory.spill': memory_spill_fraction,
            'distributed.worker.memory.pause': memory_pause_fraction,
            'distributed.worker.memory.terminate': memory_terminate_fraction,
            'distributed.comm.timeouts.tcp': tcp_timeout,
            'distributed.client.heartbeat': client_heartbeat,
            'distributed.worker.daemon': False,
        })
        
        cluster = LocalCluster(
            n_workers=n_workers,
            threads_per_worker=threads_per_worker,
            memory_limit=memory_per_worker,
            processes=processes,  # Use processes for better memory isolation
            scheduler_port=scheduler_port,
            dashboard_address=dashboard_address,
            silence_logs=True,  # Suppress Dask cluster logs
        )
        client = Client(cluster)
        
        logger.info(f"Dask dashboard: {client.dashboard_link}")
        
        return client
        
    except ImportError:
        logger.warning("Dask not available, falling back to sequential processing")
        return None