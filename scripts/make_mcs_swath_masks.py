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
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
# from src.zarr_tools import stream_process_to_zarr, initialize_zarr_store, setup_dask_client
from src.zarr_tools import setup_dask_client, initialize_zarr_store, append_chunk_to_zarr
from pyflextrkr.ft_utilities import load_config

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
    Process a time chunk of dataset to create MCS and CCS swath masks.
    
    Parameters:
    -----------
    _ds : xarray.Dataset
        Input dataset chunk with dimensions (time, cell) containing 'mcs_mask' and 'ccs_mask'.
    verbose : bool
        If True, print progress information.
    
    Returns:
    --------
    out_ds : xarray.Dataset
        Output dataset with dimensions (time, cell) containing 'mcs_mask' and 'ccs_mask' swath masks.
    """
    if verbose:
        print(f"Processing time chunk with {len(_ds.time)} time steps...")
    
    # Extract track number arrays
    mcs_mask = _ds['mcs_mask'].values  # shape (time, cell)
    ccs_mask = _ds['ccs_mask'].values  # shape (time, cell)
    
    # Create swaths and coverage for MCS
    mcs_swaths_dict, mcs_coverage_dict = create_track_swaths_and_coverage(mcs_mask)
    combined_mcs_swath = combine_swaths_with_priority(mcs_swaths_dict, mcs_coverage_dict)
    
    # Create swaths and coverage for CCS
    ccs_swaths_dict, ccs_coverage_dict = create_track_swaths_and_coverage(ccs_mask)
    combined_ccs_swath = combine_swaths_with_priority(ccs_swaths_dict, ccs_coverage_dict)
    
    if verbose:
        print(f"  ✅ Completed processing for this time chunk")

    return {
        'mcs_mask': combined_mcs_swath,
        'ccs_mask': combined_ccs_swath,
    }

def process_timechunk_wrapper_zarr(time_val, zarr_path, verbose=False):
    """
    Wrapper function for processing a time chunk by reading from zarr file.
    
    This approach avoids serialization issues by having each worker read the zarr file
    directly instead of serializing large xarray Datasets.
    
    Args:
        time_val: Time coordinate value(s) to process
        zarr_path: Path to the zarr file containing the data
        verbose: Whether to print verbose output
        
    Returns:
        tuple: (time_str, results_dict) or (time_str, None) if error
    """
    try:
        
        # Each worker opens the zarr file independently
        ds = xr.open_dataset(zarr_path, engine='zarr')
        
        # Select the specific times
        _ds = ds.sel(time=time_val)
        # Select the first time as output time value
        out_time_val = time_val[0]
        
        # Load the data into memory
        _ds = _ds.load()
        
        # Close the full dataset to free memory
        ds.close()
        
        # Process this time step
        timestep_results = process_timechunk_swath(_ds, verbose=verbose)
        
        return str(out_time_val), timestep_results
        
    except Exception as e:
        print(f"Error processing time step {time_val}: {e}")
        traceback.print_exc()
        return str(time_val), None

#--------------------------------------------------------------------------------------------------
def stream_process_to_zarr(ds, time_coords, mask_variables, output_path, template_coords, attrs,
                          client=None, logger=None, parallel=True, chunk_size_time=6, input_zarr_path=None):
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Process and write time steps in chunks
    total_processed = 0
    total_chunks = (len(time_coords) + chunk_size_time - 1) // chunk_size_time

    logger.info(f"Processing {len(time_coords)} time steps in {total_chunks} chunks of {chunk_size_time}")
    
    for chunk_idx in range(total_chunks):
        start_idx = chunk_idx * chunk_size_time
        end_idx = min((chunk_idx + 1) * chunk_size_time, len(time_coords))
        chunk_times = time_coords[start_idx:end_idx]
        
        logger.info(f"Processing chunk {chunk_idx + 1}/{total_chunks}: time steps {start_idx}-{end_idx-1}")

        # Process this chunk of time steps
        chunk_results = {}
        
        futures = []
        # Process this chunk
        future = process_timechunk_wrapper_zarr(chunk_times, input_zarr_path, verbose=False)
        futures.append(future)

        # Collect results
        for future in futures:
            time_str, result = future
            if result is not None:
                chunk_results[time_str] = result
            else:
                logger.warning(f"Skipping time step {time_str} due to processing error")

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

        import pdb; pdb.set_trace()
    return

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
    # parser.add_argument('--source', type=str, default='scream',
    #                    help='Source name (default: scream)')
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
        mask_variables = ['mcs_mask', 'ccs_mask']

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

            # Stream process with chunked zarr writing
            successful_times = stream_process_to_zarr(
                ds=ds,
                time_coords=time_coords,
                mask_variables=mask_variables,
                output_path=out_zarr,
                template_coords=ds.coords,
                attrs=attrs,
                client=client,
                logger=logger,
                parallel=parallel,
                chunk_size_time=chunk_size_time,
                input_zarr_path=in_zarr  # Pass the input zarr path for workers
            )
            import pdb; pdb.set_trace()

        except Exception as e:
            logger.error(f"Error writing chunked zarr: {e}")
            print(f"  ❌ Error writing zarr: {e}")
            return

        import pdb; pdb.set_trace()

    finally:
        # Always cleanup client
        if client and parallel:
            logger.info("Shutting down Dask client")
            client.close()

    import pdb; pdb.set_trace()

if __name__ == "__main__":
    main()