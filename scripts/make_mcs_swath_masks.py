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
                          input_zarr_path=None, batch_size=100, chunks_to_process=None):
    """
    Stream process time chunks and write results to zarr with optional parallel processing.
    
    Uses a batched submission approach to avoid overwhelming the Dask scheduler with
    too many tasks at once. Instead of submitting all chunks at once, submits them
    in batches (super-chunks), waits for each batch to complete, then moves to the next.
    
    Parameters:
    -----------
    time_coords : array-like
        Array of time coordinates to process
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
        Number of input time steps per output chunk (e.g., 6 hourly → 1 swath)
    input_zarr_path : str
        Path to input zarr file for workers to read from
    batch_size : int
        Number of chunks to submit per batch (default: 100)
    chunks_to_process : list of int, optional
        If provided, only process these chunk indices (for resume mode)
        
    Returns:
    --------
    int : Number of successfully processed chunks
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Process and write time steps in chunks
    total_processed = 0
    total_chunks = (len(time_coords) + chunk_size_time - 1) // chunk_size_time
    
    # Determine which chunks to process
    if chunks_to_process is not None:
        logger.info(f"RESUME MODE: Processing {len(chunks_to_process)} missing chunks out of {total_chunks} total")
        chunk_indices = chunks_to_process
    else:
        logger.info(f"Processing all {total_chunks} chunks")
        chunk_indices = list(range(total_chunks))

    logger.info(f"Processing {len(time_coords)} time steps in {len(chunk_indices)} chunks of {chunk_size_time}")
    logger.info(f"Each chunk will aggregate {chunk_size_time} hourly time steps into 1 swath mask")
    
    if parallel and client is not None:
        # PARALLEL MODE: Submit chunks in batches to avoid overwhelming scheduler
        total_batches = (len(chunk_indices) + batch_size - 1) // batch_size
        logger.info(f"Using batched submission: {total_batches} batches of up to {batch_size} chunks each")
        
        # Prepare chunk metadata for chunks we're actually processing
        chunk_metadata = []
        for chunk_idx in chunk_indices:
            start_idx = chunk_idx * chunk_size_time
            end_idx = min((chunk_idx + 1) * chunk_size_time, len(time_coords))
            chunk_metadata.append({
                'chunk_idx': chunk_idx,
                'start_idx': start_idx,
                'end_idx': end_idx
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
                        # Get the actual time values for this chunk for zarr writing
                        chunk_times = time_coords[start_idx:end_idx]
                        
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
            start_idx = chunk_idx * chunk_size_time
            end_idx = min((chunk_idx + 1) * chunk_size_time, len(time_coords))
            chunk_times = time_coords[start_idx:end_idx]
            
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
def check_missing_chunks(output_zarr, time_coords, chunk_size_time, mask_variables, logger=None):
    """
    Check which chunks are missing or incomplete in the output zarr file.
    
    Parameters:
    -----------
    output_zarr : str
        Path to output zarr file
    time_coords : array-like
        Expected time coordinates
    chunk_size_time : int
        Number of input time steps per output chunk
    mask_variables : list
        List of mask variable names to check
    logger : logging.Logger, optional
        Logger instance
        
    Returns:
    --------
    missing_chunks : list
        List of chunk indices that are missing or have all zeros
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    if not os.path.exists(output_zarr):
        logger.warning(f"Output file does not exist: {output_zarr}")
        return list(range((len(time_coords) + chunk_size_time - 1) // chunk_size_time))
    
    try:
        # Open zarr store
        ds_out = xr.open_zarr(output_zarr)
        
        total_chunks = (len(time_coords) + chunk_size_time - 1) // chunk_size_time
        missing_chunks = []
        
        logger.info(f"Checking {total_chunks} chunks for missing/incomplete data...")
        
        for chunk_idx in range(total_chunks):
            start_idx = chunk_idx * chunk_size_time
            end_idx = min((chunk_idx + 1) * chunk_size_time, len(time_coords))
            
            # For output, we only have 1 time per chunk (aggregated swath)
            output_time_idx = chunk_idx
            
            if output_time_idx >= len(ds_out.time):
                logger.warning(f"Chunk {chunk_idx + 1} is beyond dataset size")
                missing_chunks.append(chunk_idx)
                continue
            
            # Check if all mask variables have non-zero data
            is_missing = False
            for var_name in mask_variables:
                data_slice = ds_out[var_name].isel(time=output_time_idx).values
                if np.all(data_slice == 0):
                    is_missing = True
                    break
            
            if is_missing:
                missing_chunks.append(chunk_idx)
        
        ds_out.close()
        
        if len(missing_chunks) > 0:
            logger.warning(f"Found {len(missing_chunks)} missing/incomplete chunks: {missing_chunks[:10]}{'...' if len(missing_chunks) > 10 else ''}")
        else:
            logger.info(f"All {total_chunks} chunks are complete!")
        
        return missing_chunks
        
    except Exception as e:
        logger.error(f"Error checking zarr file: {e}")
        traceback.print_exc()
        return []

#--------------------------------------------------------------------------------------------------
def main():
    """Main function to run the make MCS swath process"""

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process co-occurrence feature overlaps')
    parser.add_argument("-c", "--config", help="yaml config file for tracking", required=True)
    parser.add_argument('--parallel', action='store_true', default=True,
                       help='Use parallel processing with Dask (default: True)')
    parser.add_argument('--no-parallel', action='store_false', dest='parallel',
                       help='Disable parallel processing')
    parser.add_argument('--workers', type=int, default=32,
                       help='Number of Dask workers (default: 32)')
    parser.add_argument('--threads-per-worker', type=int, default=1,
                       help='Number of threads per worker (default: 1)')
    parser.add_argument('--batch-size', type=int, default=100,
                       help='Number of chunks to submit per batch to avoid overwhelming scheduler (default: 100)')
    parser.add_argument('--check-missing', action='store_true',
                       help='Check for missing chunks in existing output file and exit')
    parser.add_argument('--resume', action='store_true',
                       help='Resume processing by only processing missing chunks')
    parser.add_argument('--test-steps', type=int, default=None,
                       help='Number of time steps to process for testing (default: all)')
    
    args = parser.parse_args()

    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)

    start_time = time.time()
    logger.info("Starting tracking MCS mask swath...")
    
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
    
    # Handle check-missing mode
    if args.check_missing:
        print(f"\nCHECK MISSING MODE: Analyzing existing output file...")
        try:
            ds = xr.open_zarr(in_zarr, consolidated=True, mask_and_scale=True)
            time_coords = ds["time"].values
            mask_variables = ['mcs_mask']
            chunk_size_time = 6
            
            missing_chunks = check_missing_chunks(
                output_zarr=out_zarr,
                time_coords=time_coords,
                chunk_size_time=chunk_size_time,
                mask_variables=mask_variables,
                logger=logger
            )
            
            if len(missing_chunks) > 0:
                print(f"\n⚠️  Found {len(missing_chunks)} missing chunks out of {(len(time_coords) + chunk_size_time - 1) // chunk_size_time} total")
                print(f"Missing chunk indices: {missing_chunks}")
                print(f"\nTo resume processing, run:")
                print(f"python {os.path.basename(__file__)} -c {args.config} --resume --workers {n_workers}")
            else:
                print(f"\n✅ All chunks are complete!")
            
            ds.close()
            return
            
        except Exception as e:
            print(f"❌ Error checking missing chunks: {e}")
            traceback.print_exc()
            return
    
    print(f"Parallel processing: {parallel}")
    if parallel:
        print(f"Workers: {n_workers}, Threads per worker: {threads_per_worker}")
        print(f"Batch size: {batch_size} chunks per batch")
    if args.resume:
        print(f"Resume mode: Will only process missing chunks")

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

        # Define all output variables
        mask_variables = ['mcs_mask']

        # Add processing metadata
        attrs = ds.attrs.copy()
        attrs.update({
            'processing_info': 'MCS swath mask creation',
            'processing_date': str(np.datetime64('today')),
            'processing_script': os.path.basename(__file__),
        })

        # Stream processing and writing to zarr
        print(f"\nStreaming processing and writing {len(time_coords)} time steps to zarr...")
        
        # Use a simple default chunk size for time dimension
        # With zarr-path approach, serialization is minimal regardless of chunk size
        # Chunk size only affects processing efficiency and zarr I/O
        chunk_size_time = 6
        print(f"Using default chunk_size_time={chunk_size_time} for optimal processing and zarr I/O")

        try:
            # Initialize the zarr store structure (once only)
            logger.info("Initializing zarr store...")
            initialize_zarr_store(
                output_path=out_zarr,
                time_coords=time_coords,
                mask_variables=mask_variables,
                template_coords=ds.coords,
                attrs=attrs,
                chunk_size_time=chunk_size_time
            )

            # Determine which chunks to process
            chunks_to_process = None
            if args.resume:
                logger.info("Checking for missing chunks...")
                missing_chunks = check_missing_chunks(
                    output_zarr=out_zarr,
                    time_coords=time_coords,
                    chunk_size_time=chunk_size_time,
                    mask_variables=mask_variables,
                    logger=logger
                )
                if len(missing_chunks) > 0:
                    chunks_to_process = missing_chunks
                    logger.info(f"Will process {len(missing_chunks)} missing chunks")
                else:
                    logger.info("No missing chunks found. All data is complete!")
                    return

            # Stream process with chunked zarr writing
            total_processed = stream_process_to_zarr(
                time_coords=time_coords,
                mask_variables=mask_variables,
                output_path=out_zarr,
                client=client,
                logger=logger,
                parallel=parallel,
                chunk_size_time=chunk_size_time,
                input_zarr_path=in_zarr,  # Pass the input zarr path for workers
                batch_size=batch_size,  # Number of chunks per batch
                chunks_to_process=chunks_to_process  # Only process these chunks if resume mode
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
        logger.info(f"Total processing time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
        print(f"\n✅ Processing complete!")
        print(f"   Output: {out_zarr}")
        print(f"   Total time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
        print(f"   Chunks processed: {total_processed}")


    finally:
        # Always cleanup client
        if client and parallel:
            logger.info("Shutting down Dask client")
            client.close()

if __name__ == "__main__":
    main()