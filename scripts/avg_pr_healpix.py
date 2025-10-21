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
from zarr_tools import setup_dask_client, write_zarr, zoom_level_from_nside

"""
NICAM Precipitation Temporal Resampling and Alignment Script

This script performs temporal resampling of NICAM precipitation data to align with
specific target hours (e.g., 00:00, 06:00, 12:00, 18:00 UTC).

Workflow:
1. Set up Dask distributed client for parallel processing
2. Lazily load NICAM data from intake catalog (Dask-backed)
3. Perform temporal resampling to target frequency (lazy operation)
   - Uses xarray.resample() which works efficiently with Dask
   - Aligns time bins to specified target hours
   - All operations remain lazy (no computation yet)
4. Write to Zarr format (triggers parallel computation)
   - Dask workers process chunks in parallel
   - Optimized chunking for time and spatial dimensions
   - Progress tracking via Dask dashboard

Key advantages of this approach:
- Lazy evaluation: operations build a task graph without computing
- Parallel execution: Dask distributes work across multiple workers
- Memory efficient: processes data in chunks, never loads full dataset
- Automatic alignment: pandas/xarray resampling handles time binning
"""


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


def load_dataset(catalog_path, location, source, catalog_params, vars_to_include):
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
    catalog_params : dict
        Parameters to pass to the catalog source (e.g., zoom, time)
    
    Returns:
    --------
    xarray.Dataset
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Loading dataset: {source} from {location}")
    logger.info(f"  Catalog parameters: {catalog_params}")

    # Load dataset from catalog
    cat = intake.open_catalog(catalog_path)[location]
    ds = cat[source](**catalog_params).to_dask()
    ds = ds.pipe(egh.attach_coords)

    # Select only the variables of interest
    ds = ds[vars_to_include]

    logger.info(f"  Dataset loaded: {len(ds.time)} timesteps, {len(ds.data_vars)} variables")
    logger.info(f"  Variables: {list(ds.data_vars.keys())}")
    
    return ds

def coarsen_data(ds, output_freq='6h', target_hours=[0, 6, 12, 18], rechunk_after=False):
    """
    Coarsen data temporally using time-based resampling and align to specific hours.
    
    This function performs lazy operations that work with Dask-backed datasets.
    The actual computation happens when the dataset is written to disk.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Input dataset with time dimension (can be Dask-backed)
    output_freq : str
        Target frequency for output (e.g., '6h', '3h', '12h')
    target_hours : list
        Hours of the day to align the output times to (e.g., [0, 6, 12, 18])
    rechunk_after : bool
        Whether to rechunk the dataset after resampling to optimize chunk sizes
        (default: True, recommended for better Dask performance)
    
    Returns:
    --------
    xarray.Dataset with resampled and aligned time coordinates (Dask-backed if input was)
    """
    logger = logging.getLogger(__name__)
    logger.info(f"🕒 Applying temporal resampling to {output_freq}")
    logger.info(f"   Target hours for alignment: {target_hours}")
    logger.info(f"   Original time dimension: {ds.sizes['time']} timesteps")
    
    # Get time information (only loads coordinate metadata, not actual data)
    times = pd.DatetimeIndex(ds.time.values)
    logger.info(f"   First original time: {times[0]}")
    logger.info(f"   Last original time: {times[-1]}")
    logger.info(f"   Original time spacing: {times[1] - times[0]}")
    
    # Log chunking information if dataset is chunked
    if hasattr(ds, 'chunks') and ds.chunks:
        logger.info(f"   Input chunking: {dict(ds.chunks)}")
    
    # Determine the offset to align with target hours
    # For 6h output with targets [0,6,12,18], we want bins: [0-6), [6-12), [12-18), [18-24)
    # Use 'offset' parameter to shift the bin edges
    base_hour = target_hours[0]
    offset_str = f'{base_hour}h'
    
    logger.info(f"   Using base offset: {offset_str} to align bins with target hours")
    
    # Use xarray's resample method which is built on pandas
    # This operation is LAZY with Dask - no computation happens yet!
    # The 'label' parameter controls which side of the bin gets the label
    # 'left' means the left edge (start) of each bin gets the label
    # 'closed' parameter controls which side of the interval is closed
    # 'left' means [left, right) - includes left, excludes right
    logger.info(f"   Resampling (lazy operation): freq={output_freq}, offset={offset_str}, label='left', closed='left'")
    
    ds_resampled = ds.resample(
        time=output_freq,
        offset=offset_str,
        label='left',      # Label with the left (start) edge of each interval
        closed='left'      # Left-closed intervals: [start, end)
    ).mean()
    
    logger.info(f"   After resampling: {ds_resampled.sizes['time']} timesteps")
    
    # # Rechunk if requested to optimize for subsequent operations
    # if rechunk_after:
    #     # Get spatial dimension name
    #     spatial_dim = 'cell' if 'cell' in ds_resampled.dims else 'ncol' if 'ncol' in ds_resampled.dims else None
        
    #     if spatial_dim:
    #         # Preserve spatial chunking from input, but allow time chunks to adjust
    #         original_spatial_chunks = ds.chunks.get(spatial_dim, (ds.sizes[spatial_dim],))[0]
            
    #         # For time, use reasonable chunk size based on output frequency
    #         # Estimate: if output is 6h and we want ~1 day per chunk
    #         chunk_dict = {
    #             'time': -1,  # Auto chunk time (will be set by write_zarr anyway)
    #             spatial_dim: original_spatial_chunks
    #         }
    #         logger.info(f"   Rechunking to optimize: {spatial_dim}={original_spatial_chunks}, time=auto")
    #         ds_resampled = ds_resampled.chunk(chunk_dict)
    
    # Get sample times for logging (this only loads coordinate metadata, not data)
    sample_times = pd.DatetimeIndex(ds_resampled.time.values)
    logger.info(f"   First few times: {sample_times[:min(6, len(sample_times))].strftime('%Y-%m-%d %H:%M').tolist()}")
    logger.info(f"   Last few times: {sample_times[-min(4, len(sample_times)):].strftime('%Y-%m-%d %H:%M').tolist()}")
    
    # Verify alignment
    output_hours = sample_times.hour.unique()
    output_minutes = sample_times.minute.unique()
    logger.info(f"   Unique hours in output: {sorted(output_hours)}")
    logger.info(f"   Unique minutes in output: {sorted(output_minutes)}")
    
    # Check if output is properly aligned
    misaligned = [h for h in output_hours if h not in target_hours]
    if misaligned:
        logger.warning(f"   ⚠️  Some output hours not in target list: {misaligned}")
    else:
        logger.info(f"   ✓ All output times aligned to target hours!")
    
    logger.info(f"   ✓ Resampling complete (lazy - no data computed yet)")
    
    return ds_resampled

def main():
    """Main execution function."""
    logger = setup_logging()
    
    logger.info("="*80)
    logger.info("NICAM Time Coordinate Alignment Script")
    logger.info("="*80)
    
    # Configuration
    # catalog_path = "/global/homes/f/feng045/program/hackathon/catalog/NERSC/main.yaml"
    catalog_path = "https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml"
    # location = "NERSC"
    # source = "nicam_gl11"
    # zoom = 8
    # time = "PT3H"

    location = "online"
    source = "um_glm_n2560_RAL3p3"
    catalog_params = {
        'zoom': 8,
        'time': 'PT1H',
    }
    
    # Variables to include in unified dataset
    vars_to_include = ['pr']
    vars_to_include += ['prs']  # Add more variables if needed

    # Temporal aggregation settings
    output_freq = '6h'  # Target frequency: 6-hourly
    target_hours = [0, 6, 12, 18]  # Align to these hours
    
    # Output configuration
    zoom = catalog_params['zoom']
    # output_zarr = f"/pscratch/sd/w/wcmca1/hackathon/healpix/nicam_gl11/shifted/NICAM_pr6h_z{zoom}.zarr"
    output_zarr = f"/pscratch/sd/w/wcmca1/hackathon/healpix/um_glm_n2560_RAL3p3/um_glm_n2560_RAL3p3_pr6h_z{zoom}.zarr"
    
    # Dask configuration
    use_parallel = True
    n_workers = 8
    threads_per_worker = 4
    
    # Chunk configuration for output
    chunksize_time = 28
    # chunksize_cell = None  # Let it use default
    
    logger.info("\nConfiguration:")
    logger.info(f"  Catalog: {catalog_path}")
    logger.info(f"  Source: {source}")
    logger.info(f"  Variables: {vars_to_include}")
    logger.info(f"  Output: {output_zarr}")
    logger.info(f"  Parallel: {use_parallel}")
    if use_parallel:
        logger.info(f"    Workers: {n_workers}")
        logger.info(f"    Threads per worker: {threads_per_worker}")
    
    # Step 1: Set up Dask client (before loading data for optimal parallelism)
    logger.info("\n" + "="*80)
    logger.info("Step 1: Setting up Dask client")
    logger.info("="*80)
    client = setup_dask_client(
        parallel=use_parallel,
        n_workers=n_workers,
        threads_per_worker=threads_per_worker,
        logger=logger
    )

    # Step 2: Load dataset (lazy loading with Dask)
    logger.info("\n" + "="*80)
    logger.info("Step 2: Loading dataset")
    logger.info("="*80)
    ds = load_dataset(catalog_path, location, source, catalog_params, vars_to_include)
    logger.info(f"   Dataset loaded lazily (Dask-backed): {type(ds.pr.data)}")

    # Get the HEALPix zoom level to calculate proper cell chunk size
    zoom_level = zoom_level_from_nside(ds.crs.attrs['healpix_nside'])
    chunksize_cell = 12 * 4**zoom_level

    # Step 3: Resample dataset to target frequency and alignment (lazy operation)
    logger.info("\n" + "="*80)
    logger.info("Step 3: Resampling dataset")
    logger.info("="*80)
    ds_out = coarsen_data(ds, output_freq=output_freq, target_hours=target_hours)
    logger.info(f"   Output dataset type: {type(ds_out.pr.data)}")
    
    # Step 4: Write to Zarr (this triggers parallel computation)
    logger.info("\n" + "="*80)
    logger.info("Step 4: Writing resampled dataset to Zarr")
    logger.info("="*80)
    logger.info(f"Output path: {output_zarr}")
   
    try:
        write_zarr(
            ds_out,
            output_zarr,
            client=client,
            logger=logger,
            chunksize_time=chunksize_time,
            chunksize_cell=chunksize_cell
        )
        logger.info("✓ Successfully wrote dataset to Zarr!")
        
    except Exception as e:
        logger.error(f"Error writing Zarr file: {e}")
        raise
    
    finally:
        # Clean up dask client
        if client:
            logger.info("Closing Dask client...")
            client.close()

if __name__ == "__main__":
    main()
