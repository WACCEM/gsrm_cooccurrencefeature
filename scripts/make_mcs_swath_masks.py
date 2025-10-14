import xarray as xr
import numpy as np
import pandas as pd
import cftime
import yaml
import os, glob
import time
import argparse
import logging
import traceback
import sys
import gc
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
# from src.zarr_tools import stream_process_to_zarr, initialize_zarr_store, setup_dask_client
from src.zarr_tools import setup_dask_client, initialize_zarr_store, append_chunk_to_zarr
from pyflextrkr.ft_utilities import load_config

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
def process_timechunk_swath(_ds, verbose=False):
    """
    Process a time chunk of dataset to create MCS swath masks.
    
    Parameters:
    -----------
    _ds : xarray.Dataset
        Input dataset chunk with dimensions (time, cell) containing 'mcs_mask'.
    verbose : bool
        If True, print progress information.
    
    Returns:
    --------
    out_ds : xarray.Dataset
        Output dataset with dimensions (time, cell) containing 'mcs_mask' swath mask.
    """
    if verbose:
        print(f"Processing time chunk with {len(_ds.time)} time steps...")
    
    # Extract track number arrays
    mcs_mask = _ds['mcs_mask'].values  # shape (time, cell)
    
    # Create swaths and coverage for MCS
    mcs_swaths_dict, mcs_coverage_dict = create_track_swaths_and_coverage(mcs_mask)
    combined_mcs_swath = combine_swaths_with_priority(mcs_swaths_dict, mcs_coverage_dict)
    
    if verbose:
        print(f"  ✅ Completed processing for this time chunk")

    return {
        'mcs_mask': combined_mcs_swath,
    }

def process_timechunk_wrapper_zarr(start_idx, end_idx, zarr_path, verbose=False):
    """
    Wrapper function for processing a time chunk by reading from zarr file.
    
    This approach avoids serialization issues by:
    1. Having each worker read the zarr file directly
    2. Using integer indices instead of passing time arrays (reduces graph size)
    
    Args:
        start_idx: Start index in the time dimension
        end_idx: End index in the time dimension (exclusive)
        zarr_path: Path to the zarr file containing the data
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
        
        # Close the full dataset to free memory
        ds.close()
        
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
                          client=None, logger=None, parallel=True, chunk_size_time=6, 
                          input_zarr_path=None, batch_size=100,
                          time_groups=None, output_time_coords=None):
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
    client : dask.distributed.Client, optional
        Dask client for parallel processing
    logger : logging.Logger, optional
        Logger instance
    parallel : bool
        Whether to use parallel processing
    chunk_size_time : int
        Aggregation window - number of input hourly time steps to aggregate into 1 swath mask (e.g., 6)
    input_zarr_path : str
        Path to input zarr file for workers to read from
    batch_size : int
        Number of output chunks to submit per batch (default: 100)
    time_groups : dict, optional
        Dictionary mapping aligned output times to lists of input time indices
    output_time_coords : array-like, optional
        Array of aligned output time coordinates
        
    Returns:
    --------
    int : Number of successfully processed chunks
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Use time groups if provided, otherwise fall back to simple chunking
    if time_groups is not None and output_time_coords is not None:
        # Time-aligned processing: use time groups
        total_chunks = len(output_time_coords)
        output_times_list = list(output_time_coords)
        chunk_indices = list(range(total_chunks))
        
        logger.info(f"Processing {len(time_coords)} input time steps into {total_chunks} aligned time groups")
        logger.info(f"Each chunk aggregates multiple hourly time steps into 1 swath mask at standard hours")
        
    else:
        # Legacy simple chunking (every N time steps)
        total_chunks = (len(time_coords) + chunk_size_time - 1) // chunk_size_time
        chunk_indices = list(range(total_chunks))

        logger.info(f"Processing {len(time_coords)} time steps in {total_chunks} chunks of {chunk_size_time}")
        logger.info(f"Each chunk will aggregate {chunk_size_time} hourly time steps into 1 swath mask")
    
    # Process and write time steps in chunks
    total_processed = 0
    
    if parallel and client is not None:
        # PARALLEL MODE: Submit chunks in batches to avoid overwhelming scheduler
        total_batches = (len(chunk_indices) + batch_size - 1) // batch_size
        logger.info(f"Using batched submission: {total_batches} batches of up to {batch_size} chunks each")
        
        # Prepare chunk metadata for chunks we're actually processing
        chunk_metadata = []
        for chunk_idx in chunk_indices:
            if time_groups is not None and output_time_coords is not None:
                # Time-aligned processing: get indices from time groups
                aligned_time = output_times_list[chunk_idx]
                # Convert numpy datetime64 to pandas Timestamp for dictionary lookup
                aligned_time_pd = pd.Timestamp(aligned_time)
                time_indices = time_groups[aligned_time_pd]
                start_idx = int(time_indices[0])
                end_idx = int(time_indices[-1]) + 1
                output_time = aligned_time
            else:
                # Legacy simple chunking
                start_idx = chunk_idx * chunk_size_time
                end_idx = min((chunk_idx + 1) * chunk_size_time, len(time_coords))
                output_time = time_coords[start_idx]
            
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
            if time_groups is not None and output_time_coords is not None:
                # Time-aligned processing: get indices from time groups
                aligned_time = output_times_list[chunk_idx]
                # Convert numpy datetime64 to pandas Timestamp for dictionary lookup
                aligned_time_pd = pd.Timestamp(aligned_time)
                time_indices = time_groups[aligned_time_pd]
                start_idx = int(time_indices[0])
                end_idx = int(time_indices[-1]) + 1
                output_time_val = aligned_time
            else:
                # Legacy simple chunking
                start_idx = chunk_idx * chunk_size_time
                end_idx = min((chunk_idx + 1) * chunk_size_time, len(time_coords))
                output_time_val = time_coords[start_idx]
            
            logger.info(f"Processing chunk {chunk_idx + 1}/{total_chunks}: time steps {start_idx}-{end_idx-1}")
            
            chunk_results = {}
            time_str, result = process_timechunk_wrapper_zarr(start_idx, end_idx, input_zarr_path, verbose=False)
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
    parser.add_argument("-c", "--config", help="yaml config file for tracking", required=True)
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
    
    # Define output variables (default only do mcs_mask)
    mask_variables = ['mcs_mask']
    
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
    
    # Get source name from root path (e.g., /pscratch/sd/w/wcmca1/hackathon/mcs/scream/)
    source_name = os.path.basename(os.path.normpath(root_path))

    # Output paths
    out_dir = "/pscratch/sd/w/wcmca1/hackathon/mcs_masks/"
    out_basename = f"{source_name}_mcs_masks_hp{zoom}.zarr"
    out_zarr = f"{out_dir}{out_basename}"
    os.makedirs(out_dir, exist_ok=True)

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
                chunk_size_time=zarr_chunk_size_time  # Zarr storage chunks (28 = 1 week)
            )

            # Stream process with chunked zarr writing
            total_processed = stream_process_to_zarr(
                time_coords=time_coords,  # Pass input time coords for processing
                mask_variables=mask_variables,
                output_path=out_zarr,
                client=client,
                logger=logger,
                parallel=parallel,
                chunk_size_time=aggregation_window,  # Aggregation window (e.g., 6 hours)
                input_zarr_path=in_zarr,  # Pass the input zarr path for workers
                batch_size=batch_size,  # Number of chunks per batch
                time_groups=time_groups,  # Mapping of aligned times to input indices
                output_time_coords=output_time_coords  # Aligned output time coordinates
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