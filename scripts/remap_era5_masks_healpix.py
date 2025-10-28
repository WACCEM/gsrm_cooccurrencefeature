import xarray as xr
import numpy as np
import os
import glob
import time
import logging
import sys
from pathlib import Path
import healpy as hp
from functools import partial
import zarr

# Add src directory to path for zarr_tools import
sys.path.append(str(Path(__file__).parent.parent / 'src'))
from zarr_tools import setup_dask_client

def setup_logging():
    """
    Set the logging message level

    Args:
        None.

    Returns:
        None.
    """
    logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO)


def get_datasets(files_ar, files_tc, files_etc, chunk_netcdf=None, parallel=False, logger=None):
    """
    Load datasets from files using optimized multi-file reading
    
    Args:
        files_ar: list
            List of AR NetCDF files
        files_tc: list
            List of TC NetCDF files
        files_etc: list
            List of ETC NetCDF files
        chunk_netcdf: dict, optional
            Dictionary specifying chunk sizes for xarray.open_mfdataset
        parallel: bool
            Whether to use parallel processing
        logger: logging.Logger
            Logger for status messages
            
    Returns:
        tuple: (ds_ar, ds_tc, ds_etc) datasets
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    def process_dataset(files, name):
        """Helper to process each dataset with optimized multi-file reading"""
        logger.info(f"Reading {name} files using open_mfdataset...")
        
        # Use open_mfdataset for efficient multi-file reading
        # Use chunking to avoid loading entire dataset into memory        
        try:
            ds = xr.open_mfdataset(
                files,
                combine='nested',           # Faster than 'by_coords' when time order is known
                concat_dim='time',          # Concatenate along time
                parallel=parallel,          # Use dask for parallel file reading
                engine='netcdf4',           # Explicit engine for consistency
                mask_and_scale=False,       # Skip scaling for performance
                decode_times=True,          # Keep time decoding for proper sorting
                chunks=chunk_netcdf,          # CRITICAL: Enable chunking
            )
            
            # Check for duplicate time values (less common with open_mfdataset)
            time_values = ds['time'].values
            unique_times, unique_indices = np.unique(time_values, return_index=True)
            
            if len(unique_times) < len(time_values):
                logger.warning(f"Found {len(time_values) - len(unique_times)} duplicate time values in {name} files")
                # Remove duplicates by keeping first occurrence
                ds = ds.isel(time=sorted(unique_indices))
            
            logger.info(f"Finished reading {name} files: {len(ds.time)} timesteps")
            logger.info(f"  Initial dataset chunks: {dict(ds.chunks)}")
            
            # Rechunk immediately after loading
            # This is more efficient than rechunking after remapping because:
            # 1. Smaller task graph during remapping (for larger time chunks)
            # 2. No expensive rechunk operation before zarr write
            # 3. Better memory locality during computation
            if chunk_netcdf is not None:
                logger.info(f"Rechunking {name} dataset to target chunks: {chunk_netcdf}")
                ds = ds.chunk(chunk_netcdf)
                logger.info(f"  Rechunked dataset chunks: {dict(ds.chunks)}")
            
        except Exception as e:
            logger.error(f"Error opening {name} files with open_mfdataset: {e}")
            raise
    
        return ds
    
    # Process each dataset using our helper function
    ds_ar = process_dataset(files_ar, "AR")
    ds_tc = process_dataset(files_tc, "TC")
    ds_etc = process_dataset(files_etc, "ETC")
    
    return ds_ar, ds_tc, ds_etc

def combine_masks(ds_ar, ds_tc, ds_etc, logger=None):
    """
    Combine AR, TC, ETC tracking datasets.
    
    Args:
        ds_ar: xarray.Dataset
            AR tracking dataset
        ds_tc: xarray.Dataset
            TC tracking dataset
        ds_etc: xarray.Dataset
            ETC tracking dataset
        logger: logging.Logger, optional
            Logger for status messages
            
    Returns:
        xarray.Dataset: Combined dataset with all masks
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    drop_var_list = ['ETC_binary_tag']
    rename_dict = {
        # 'AR_binary_tag': 'ar_mask',
        # 'TC_binary_tag': 'tc_mask',
        # 'ETC_binary_tag': 'etc_mask',
        'AR_count_index': 'ar_mask',
        'TC_count_index': 'tc_mask',
        'ETC_int_tag': 'etc_mask',
        'longitude': 'lon',
        'latitude': 'lat',
    }
    
    # Find common time range more efficiently
    # Convert to sets for faster intersection
    logger.info("Finding common time range across datasets...")
    common_times = sorted(set(ds_ar['time'].values)
                         .intersection(set(ds_tc['time'].values))
                         .intersection(set(ds_etc['time'].values)))
    if not common_times:
        logger.warning("No common time values between all datasets!")
        return None
    else:
        logger.info(f"Found {len(common_times)} common time points")
        # Select only the common times in all datasets
        ds_ar = ds_ar.sel(time=common_times)
        ds_tc = ds_tc.sel(time=common_times)
        ds_etc = ds_etc.sel(time=common_times)

    # Merge the datasets
    logger.info("Merging datasets...")
    ds = xr.merge([ds_ar, ds_tc, ds_etc], combine_attrs='drop_conflicts', compat='override')
    logger.info(f"Successfully merged datasets with {len(common_times)} common time points")

    # Rename variables, drop unwanted ones in the DataSet
    ds = ds.rename(rename_dict).drop_vars(drop_var_list, errors='ignore')

    # TODO: Modify global attributes if needed
    # ds.attrs['history'] = f"Created on {time.ctime()} by combining tracking data"
    
    return ds


def fix_coords(ds, lat_dim="lat", lon_dim="lon", roll=False):
    """
    Fix coordinates in a dataset:
    1. Convert longitude from -180/+180 to 0-360 range (optional)
    2. Roll dataset to start at longitude 0 (optional)
    3. Ensure coordinates are in ascending order
    
    Parameters:
    -----------
    ds : xarray.Dataset or xarray.DataArray
        Dataset with lat/lon coordinates
    lat_dim : str, optional
        Name of latitude dimension, default "lat"
    lon_dim : str, optional
        Name of longitude dimension, default "lon"
    roll : bool, optional, default=False
        If True, convert longitude from -180/+180 to 0-360, and roll the dataset to start at longitude 0
        
    Returns:
    --------
    xarray.Dataset or xarray.DataArray
        Dataset with fixed coordinates
    """
    if roll:
        # Find where longitude crosses from negative to positive (approx. where lon=0)
        lon_0_index = (ds[lon_dim] < 0).sum().item()
        
        # Create indexers for the roll
        lon_indices = np.roll(np.arange(ds.sizes[lon_dim]), -lon_0_index)
        
        # Roll dataset and convert longitudes to 0-360 range
        ds = ds.isel({lon_dim: lon_indices})
        lon360 = xr.where(ds[lon_dim] < 0, ds[lon_dim] + 360, ds[lon_dim])
        ds = ds.assign_coords({lon_dim: lon360})
    
    # Ensure latitude and longitude are in ascending order if needed
    if np.all(np.diff(ds[lat_dim].values) < 0):
        ds = ds.isel({lat_dim: slice(None, None, -1)})
    if np.all(np.diff(ds[lon_dim].values) < 0):
        ds = ds.isel({lon_dim: slice(None, None, -1)})
    
    return ds


def is_valid(ds, tolerance=0.1):
    """
    Limit extrapolation distance to a certain tolerance.
    This is useful for preventing extrapolation of regional data to global HEALPix grid.

    Args:
        ds (xarray.Dataset):
            The dataset containing latitude and longitude coordinates.
        tolerance (float): default=0.1
            The maximum allowed distance in [degrees] for extrapolation.

    Returns:
        xarray.DataSet.
    """
    return (np.abs(ds.lat - ds.lat_hp) < tolerance) & (np.abs(ds.lon - ds.lon_hp) < tolerance)


def calculate_healpix_tolerance(zoom_level):
    """
    Calculate appropriate tolerance for is_valid function based on HEALPix zoom level.
    Returns approximately one grid cell size in degrees.
    
    Args:
        zoom_level (int): HEALPix zoom level
        
    Returns:
        float: Tolerance in degrees
    """
    # Calculate nside from zoom level (nside = 2^zoom)
    # nside determines HEALPix resolution - each increase in zoom doubles the resolution
    nside = 2 ** zoom_level
    
    # Calculate approximate pixel size in degrees
    # Mathematical derivation:
    # - Sphere has total area of 4π steradians (= 4π × (180/π)² sq. degrees)
    # - HEALPix divides sphere into 12 × nside² equal-area pixels
    # - Each pixel has area = 4π × (180/π)² / (12 × nside²) sq. degrees
    # - Linear size = √(pixel area) ≈ 58.6 / nside degrees
    # This gives approximately the angular width of one HEALPix cell
    pixel_size_degrees = 58.6 / nside
    
    return pixel_size_degrees

def remap_to_healpix_and_save(ds, zoom, out_zarr, 
                                chunksize_cell, chunksize_time,
                                client=None, logger=None):
    """
    Remap a dataset to HEALPix grid and save as Zarr.
    
    Args:
        ds : xarray.Dataset
            Input dataset to remap.
        zoom : int
            HEALPix zoom level.
        out_zarr : str
            Output Zarr file path.
        chunksize_cell : int
            Chunk size for cell dimension.
        chunksize_time : int
            Chunk size for time dimension.
        client : (dask.distributed.Client, optional)
            Dask client.
        logger : (logging.Logger, optional)
            Logger for debug information.

    Returns:
        xarray.Dataset: The remapped HEALPix dataset
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Generate HEALPix coordinates using healpy
    nside = 2**zoom
    npix = hp.nside2npix(nside)
    logger.info(f"Generating HEALPix coordinates for zoom={zoom} (nside={nside}, npix={npix})")
    
    # Get lon/lat for all HEALPix cells using healpy
    lon_vals, lat_vals = hp.pix2ang(
        nside=nside, ipix=np.arange(npix), lonlat=True, nest=True
    )
    
    # Convert numpy arrays to xarray DataArrays with proper coordinates
    # This is critical for performance of the remapping using xarray's .sel() method
    # The key is to have 'lon', 'lat' coordinates (without _hp) plus 'lon_hp', 'lat_hp'
    cell_coord = np.arange(npix)
    
    # Create DataArrays matching the catalog structure exactly
    lon_hp = xr.DataArray(
        lon_vals,
        dims=['cell'],
        coords={
            'cell': cell_coord,
            'lat': ('cell', lat_vals),      # Coordinate without _hp suffix
            'lon': ('cell', lon_vals),      # Coordinate without _hp suffix  
            'lon_hp': ('cell', lon_vals),   # Self-reference for tolerance check
        },
        attrs={
            'units': 'degree_east',
            'standard_name': 'longitude',
            'axis': 'X',
        }
    )
    
    lat_hp = xr.DataArray(
        lat_vals,
        dims=['cell'],
        coords={
            'cell': cell_coord,
            'lat': ('cell', lat_vals),      # Coordinate without _hp suffix
            'lon': ('cell', lon_vals),      # Coordinate without _hp suffix
            'lat_hp': ('cell', lat_vals),   # Self-reference for tolerance check
        },
        attrs={
            'units': 'degree_north',
            'standard_name': 'latitude',
            'axis': 'Y',
        }
    )
    
    logger.info(f"Created HEALPix DataArrays with coordinates for efficient indexing")
    
    # Make sure coordinates are fixed before remapping
    ds = fix_coords(ds)

    # Calculate appropriate tolerance based on zoom level
    tolerance = calculate_healpix_tolerance(zoom)
    logger.info(f"Using HEALPix tolerance of {tolerance:.4f}° at zoom level {zoom}")
    
    logger.info("Applying nearest neighbor remapping to HEALPix grid...")
    logger.info(f"Remapping dataset with {len(ds.time)} timesteps in chunks")
    
    fill_value = 0
    
    # Remap the entire dataset at once (xarray handles chunking internally)
    dsout_hp = ds.sel(
        lon=lon_hp, lat=lat_hp, method="nearest",
    ).where(partial(is_valid, tolerance=tolerance), fill_value)

    # Drop lat/lon coordinates (not needed in HEALPix)
    dsout_hp = dsout_hp.drop_vars(["lat_hp", "lon_hp", "lat", "lon"])
    
    # Update global attributes
    dsout_hp.attrs['Title'] = f"HEALPix remapped tracking mask data (zoom={zoom})"
    dsout_hp.attrs['zoom'] = zoom
    dsout_hp.attrs["Created_on"] = time.ctime(time.time())
    dsout_hp.attrs['processing_script'] = os.path.basename(__file__)
    
    # Check current chunking - should already be optimal from input rechunking
    logger.info(f"HEALPix dataset dimensions: {dict(dsout_hp.sizes)}")
    logger.info(f"Current HEALPix chunks: {dict(dsout_hp.chunks)}")
    
    # Only rechunk if necessary (e.g., if chunks don't match desired output)
    current_time_chunks = dsout_hp.chunks.get('time', [])
    needs_rechunk = False
    
    # Check if time chunking needs adjustment
    if current_time_chunks:
        max_time_chunk = max(len(c) if isinstance(c, (list, tuple)) else c for c in [current_time_chunks])
        if max_time_chunk != chunksize_time:
            needs_rechunk = True
            logger.info(f"Time chunks ({max_time_chunk}) don't match target ({chunksize_time}), rechunking...")
    
    if needs_rechunk:
        logger.info(f"Rechunking HEALPix output to: time={chunksize_time}, cell={chunksize_cell}")
        chunked_hp = dsout_hp.chunk({
            "time": chunksize_time, 
            "cell": chunksize_cell, 
        })
    else:
        # Data is already optimally chunked, just ensure cell dimension is right
        logger.info(f"Dataset already has optimal time chunking, adjusting cell chunks only")
        chunked_hp = dsout_hp.chunk({"cell": chunksize_cell})
    
    # Report final chunking scheme
    n_time_chunks = len(chunked_hp.time) // chunksize_time + (1 if len(chunked_hp.time) % chunksize_time else 0)
    n_cell_chunks = npix // chunksize_cell + (1 if npix % chunksize_cell else 0)
    logger.info(f"Final chunks: {n_time_chunks} (time) × {n_cell_chunks} (cell) = {n_time_chunks * n_cell_chunks}")


    # ---------- WRITE HEALPIX ZARR OUTPUT ----------
    logger.info(f"Starting HEALPix Zarr write to: {out_zarr}")
    
    # Configure compression for all variables
    # Blosc with zstd + bitshuffle is optimal for integer mask data
    # - cname='zstd': Use Zstd compression algorithm
    # - clevel=3: Compression level (1-9, 3 is good balance)
    # - shuffle=2: BITSHUFFLE - very effective for integer arrays with repeated values
    compressor = zarr.Blosc(cname='zstd', clevel=3, shuffle=2)
    
    encoding = {}
    for var in chunked_hp.data_vars:
        encoding[var] = {'compressor': compressor}
    
    # Create a delayed task for Zarr writing
    write_task = chunked_hp.to_zarr(
        out_zarr,
        mode="w",
        consolidated=True,  # Enable for better performance when reading
        encoding=encoding,
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
        try:
            memory_usage = client.run(lambda: psutil.Process().memory_info().rss / 1e9)
            logger.info(f"Current memory usage across workers (GB): {memory_usage}")
        except Exception as e:
            logger.warning(f"Could not get memory usage: {e}")
                
        try:
            # Compute with progress tracking
            future = client.compute(write_task)
            logger.info("Writing HEALPix Zarr (this may take a while)...")
            logger.info("Monitor progress at Dask dashboard")
            progress(future)  # Shows a progress bar in notebooks or detailed progress in terminals

            result = future.result()
            logger.info("HEALPix Zarr write completed successfully")
        except Exception as e:
            logger.error(f"HEALPix Zarr write failed: {str(e)}")
            raise
        finally:
            # Restore original log level
            shuffle_logger.setLevel(original_level)
    else:
        # Compute locally if no client
        logger.info("Computing Zarr write (no distributed client)...")
        write_task.compute()

    logger.info(f"✅ HEALPix conversion complete: {out_zarr}")
    
    return dsout_hp

def main():
    """Main function to run the remap masks"""
    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)

    start_time = time.time()
    logger.info("Starting remap masks ...")

    # Define parameters
    source_name = "ERA5"
    zoom = 8
    version = "v1"
    parallel = True
    
    # Dask setup parameters
    n_workers = 8
    threads_per_worker = 4
    
    chunksize_cell = 12 * 4**zoom
    chunksize_time = 28

    # Bryce's original ERA5 tracking files
    # in_dir = "/pscratch/sd/b/beharrop/kmscale_hackathon/ERA5_tracking/"
    in_dir = "/pscratch/sd/b/beharrop/kmscale_hackathon/hackathon_pre/era5_tracking/"
    dir_ar = f"{in_dir}"
    dir_tc = f"{in_dir}"
    dir_etc = f"{in_dir}"
    basename_ar = f"AR_tracks_era5_*.nc"
    basename_tc = f"TC_tracks_era5_*.nc"
    basename_etc = f"ETC_test_tracks_era5_*.nc"

    out_dir = "/pscratch/sd/w/wcmca1/hackathon/all_masks/"
    out_basename = f"{source_name}_AR_TC_ETC_hp{zoom}_{version}.zarr"
    out_zarr = f"{out_dir}{out_basename}"
    os.makedirs(out_dir, exist_ok=True)

    # Specify chunking for netCDF reading
    chunk_netcdf = {
        'time': chunksize_time,  # Set time chunk to match output
        'latitude': -1,          # No chunking on latitude dimension
        'longitude': -1,         # No chunking on longitude dimension
    }

    # Setup Dask client
    client = setup_dask_client(parallel=parallel, n_workers=n_workers, threads_per_worker=threads_per_worker, logger=logger)

    try:
        # Find input files
        files_ar = sorted(glob.glob(f"{dir_ar}{basename_ar}"))
        files_tc = sorted(glob.glob(f"{dir_tc}{basename_tc}"))
        files_etc = sorted(glob.glob(f"{dir_etc}{basename_etc}"))
        logger.info(f"Number of AR files: {len(files_ar)}")
        logger.info(f"Number of TC files: {len(files_tc)}")
        logger.info(f"Number of ETC files: {len(files_etc)}")

        # Load datasets
        ds_ar, ds_tc, ds_etc = get_datasets(files_ar, files_tc, files_etc, chunk_netcdf=chunk_netcdf, parallel=False, logger=logger)

        # Combine datasets
        ds = combine_masks(ds_ar, ds_tc, ds_etc, logger=logger)

        # Remap to HEALPix and save
        dsout_hp = remap_to_healpix_and_save(ds, zoom, out_zarr, 
                                        chunksize_cell, chunksize_time,
                                        client=client, logger=logger)
        
        # Cleanup
        ds_ar.close()
        ds_tc.close()
        ds_etc.close()
        ds.close()
        dsout_hp.close()
        
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