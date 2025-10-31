import xarray as xr
import numpy as np
import pandas as pd
import cftime
import os, glob
import time
import logging
import sys
import yaml
import argparse
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from src.zarr_tools import zoom_level_from_nside, write_zarr, setup_dask_client
from src.utilities import convert_to_matching_calendar

#-------------------------------------------------------------------
def combine_masks(ds_mcs, ds_ar, ds_tc, ds_etc, client=None, out_zarr=None, logger=None):
    """
    Combine tracking datasets and write to Zarr store.
    
    Args:
        ds_mcs: xarray.Dataset
            MCS tracking dataset
        ds_ar: xarray.Dataset
            AR tracking dataset
        ds_tc: xarray.Dataset
            TC tracking dataset
        ds_etc: xarray.Dataset
            ETC tracking dataset
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

    drop_var_list = ['pr', 'ETC_binary_tag', 'sfcWind']
    rename_dict = {
        'AR_count_index': 'ar_mask',
        'TC_count_index': 'tc_mask',
        # 'ETC_count_index': 'etc_mask',
        'ETC_int_tag': 'etc_mask',
    }

    # Check the calendar type of the time coordinate in ds_ar
    calendar = ds_ar['time'].dt.calendar

    # Convert ds_mcs time coordinate to match ds_ar calendar
    if calendar not in ['proleptic_gregorian', 'gregorian', 'standard']:
        logger.info(f"Converting ds_mcs time from standard calendar to {calendar} calendar")
        converted_times = convert_to_matching_calendar(ds_mcs['time'].values, calendar)
        
        # Replace the time values in ds_mcs with the converted ones
        ds_mcs = ds_mcs.assign_coords(time=converted_times)
        logger.info(f"Successfully converted ds_mcs time to {calendar} calendar")
    else:
        logger.info(f"Both datasets use standard calendar: {calendar}, no conversion needed")
    
    # Find common time range across all three datasets
    common_times = sorted(set(ds_mcs['time'].values)
                         .intersection(set(ds_ar['time'].values))
                         .intersection(set(ds_tc['time'].values))
                         .intersection(set(ds_etc['time'].values)))
    if not common_times:
        logger.warning("No common time values between all datasets!")
        return None
    else:
        # Select only the common times in all datasets
        ds_mcs = ds_mcs.sel(time=common_times)
        ds_ar = ds_ar.sel(time=common_times)
        ds_tc = ds_tc.sel(time=common_times)
        ds_etc = ds_etc.sel(time=common_times)

    # Fix for lat/lon coordinates issue: ensure consistent treatment
    datasets = [ds_mcs, ds_ar, ds_tc, ds_etc]
    # Fix dimensions and coordinates for consistent merging
    for i, ds in enumerate(datasets):
        # Rename 'ncol' to 'cell' if it exists to standardize dimensions
        if 'ncol' in ds.dims:
            logger.info(f"Renaming dimension 'ncol' to 'cell' in dataset {i}")
            datasets[i] = ds.rename({'ncol': 'cell'})
        
        # Drop lat and lon variables from all datasets
        for var in ['lat', 'lon']:
            if var in ds.variables:
                logger.info(f"Dropping {var} from dataset {i}")
                datasets[i] = datasets[i].drop_vars(var)

    # Merge the datasets
    ds = xr.merge(datasets, combine_attrs='drop_conflicts', compat='override')
    logger.info(f"Successfully merged datasets with {len(common_times)} common time points")

    # Rename variables, drop unwanted ones in the DataSet
    ds = ds.rename(rename_dict).drop_vars(drop_var_list, errors='ignore')

    # Get the HEALPix zoom level to calculate proper cell chunk size
    zoom_level = zoom_level_from_nside(ds.crs.attrs['healpix_nside'])
    chunksize_cell = 12 * 4**zoom_level
    chunksize_time = 28

    # Rechunk to ensure consistent chunking across all variables
    ds = ds.chunk({'time': chunksize_time, 'cell': chunksize_cell})
    logger.info(f"After rechunk: {dict(ds.chunks)}")

    # Write to Zarr
    if out_zarr:
        write_zarr(ds, out_zarr, client=client, logger=logger)
    
    return ds

def setup_logging():
    """
    Set the logging message level

    Args:
        None.

    Returns:
        None.
    """
    logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO)

def load_config(config_file, catalog_source):
    """
    Load configuration from YAML file for a specific catalog source.
    
    Args:
        config_file: str
            Path to the YAML configuration file
        catalog_source: str
            The catalog source key to load configuration for
            
    Returns:
        dict: Configuration dictionary for the specified source
    """
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    if catalog_source not in config:
        raise ValueError(f"Catalog source '{catalog_source}' not found in config file. "
                        f"Available sources: {list(config.keys())}")
    
    return config[catalog_source]

def get_datasets(dir_mcs, files_ar, files_tc, files_etc, parallel=False, logger=None):
    """
    Load datasets from files
    
    Args:
        dir_mcs: str
            Directory containing MCS zarr store
        files_ar: list
            List of AR NetCDF files
        files_tc: list
            List of TC NetCDF files
        files_etc: list
            List of ETC NetCDF files
        parallel: bool
            Whether to use parallel processing
        logger: logging.Logger
            Logger for status messages
            
    Returns:
        tuple: (ds_mcs, ds_ar) datasets
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    def process_dataset(files, name):
        """Helper to process each dataset with consistent error handling"""
        logger.info(f"Reading {name} files...")
        
        # Process each file individually, then concatenate
        datasets = []
        for file in files:
            try:
                ds_single = xr.open_dataset(
                    file,
                    chunks={},
                    mask_and_scale=False,
                )
                datasets.append(ds_single)
            except Exception as e:
                logger.warning(f"Error opening {file}: {e}")
                continue
        
        if not datasets:
            raise ValueError(f"Could not open any {name} files")
        
        # Concatenate along time dimension
        ds = xr.concat(datasets, dim="time")
        
        # Sort time values to ensure monotonic order
        logger.info(f"Sorting {name} dataset by time")
        ds = ds.sortby('time')
        
        # Check for and remove duplicate time values
        _, index = np.unique(ds['time'].values, return_index=True)
        if len(index) < len(ds['time']):
            logger.warning(f"Found {len(ds['time']) - len(index)} duplicate time values in {name} files")
            ds = ds.isel(time=sorted(index))
    
        logger.info(f"Finished reading {name} files.")
        return ds
        
    # Process each dataset using our helper function
    ds_ar = process_dataset(files_ar, "AR")
    ds_tc = process_dataset(files_tc, "TC")
    ds_etc = process_dataset(files_etc, "ETC")
    
    # Read MCS file
    logger.info("Reading MCS file...")
    ds_mcs = xr.open_zarr(
        dir_mcs,
        consolidated=True,
        mask_and_scale=False,
    )
    logger.info(f"Finished reading MCS file")
    
    return ds_mcs, ds_ar, ds_tc, ds_etc

def main():
    """Main function to run the mask combination process"""

    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    start_time = time.time()
    logger.info("Starting tracking mask combination...")

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process combined feature masks')
    parser.add_argument("-c", "--config", help="yaml config file for processing", required=True)
    parser.add_argument('--source', type=str, help='Source name', required=True)

    # Configuration from arguments
    args = parser.parse_args()
    source = args.source
    config_file = args.config

    # Load configuration for the specified source
    config = load_config(config_file, source)
    dir_te = config.get("dir_te")
    source_name = config.get("source_name")
    source_te = config.get("source_te")
    source_res = config.get("source_res")
    
    # Configuration parameters
    zoom = 8
    version = 'v1'
    parallel = True
    n_workers = 16
    threads_per_worker = 4
    
    # Input/output paths
    # source_name = "scream"
    # source_res = f"ne120_inst_ivt_hp{zoom}"
    # dir_mcs = f"/pscratch/sd/w/wcmca1/scream-cess-healpix/mcs_tracking_hp9/mcstracking/{source_name}2D_hrly_mcsmask_hp8_v1.zarr"
    # dir_te = f"/pscratch/sd/b/beharrop/kmscale_hackathon/hackathon_pre/scream_1year_test/"
    # basename_ar = f"AR_tracks_{source_name}2D_ne120_inst_ivt_hp8."
    # basename_tc = f"TC_tracks_{source_name}2D_ne120_inst_ivt_hp8."
    # basename_etc = f"ETC_tracks_{source_name}2D_ne120_inst_ivt_hp8."

    # source_name = "icon_d3hp003"
    # source_res = f"hp{zoom}_PT6H"
    # dir_mcs = f"/pscratch/sd/w/wcmca1/hackathon/mcs/{source_name}/mcstracking/icon_hrly_mcsmask_hp{zoom}_v1.zarr"
    # dir_te = f"/pscratch/sd/b/beharrop/kmscale_hackathon/hackathon_pre/{source_name}_1year_testpy/"
    # basename_ar = f"AR_tracks_{source_name}_{source_res}."
    # basename_tc = f"TC_tracks_{source_name}_{source_res}."
    # basename_etc = f"ETC_tracks_{source_name}_{source_res}."

    # source_name = "casesm2_10km_nocumulus"
    # source_res = f"hp{zoom}_H"
    # dir_mcs = f"/pscratch/sd/w/wcmca1/hackathon/mcs/{source_name}/mcstracking/casesm2_hrly_mcsmask_hp{zoom}_v1.zarr"
    dir_mcs = f"/pscratch/sd/w/wcmca1/hackathon/mcs_masks/{source_name}_mcs_masks_hp{zoom}.zarr"
    basename_ar = f"AR_tracks_{source_te}_{source_res}."
    basename_tc = f"TC_tracks_{source_te}_{source_res}."
    basename_etc = f"ETC_test_tracks_{source_te}_{source_res}."

    # Output paths
    out_dir = "/pscratch/sd/w/wcmca1/hackathon/all_masks/"
    out_basename = f"{source_name}_allmasks_hp{zoom}_{version}.zarr"
    out_zarr = f"{out_dir}{out_basename}"
    os.makedirs(out_dir, exist_ok=True)
    
    # Setup Dask client
    client = setup_dask_client(parallel=parallel, n_workers=n_workers, threads_per_worker=threads_per_worker, logger=logger)
    
    try:
        # Find input files
        files_ar = sorted(glob.glob(f"{dir_te}{basename_ar}*.nc"))
        files_tc = sorted(glob.glob(f"{dir_te}{basename_tc}*.nc"))
        files_etc = sorted(glob.glob(f"{dir_te}{basename_etc}*.nc"))
        logger.info(f"Number of AR files: {len(files_ar)}")
        logger.info(f"Number of TC files: {len(files_tc)}")
        logger.info(f"Number of ETC files: {len(files_etc)}")

        # Check MCS directory
        if not os.path.exists(dir_mcs):
            logger.error(f"MCS directory does not exist: {dir_mcs}")
            sys.exit(1)
        # Check TE files
        if len(files_ar) == 0:
            logger.error(f"Missing AR files")
            sys.exit(1)
        if len(files_tc) == 0:
            logger.error(f"Missing TC files")
            sys.exit(1)
        if len(files_etc) == 0:
            logger.error(f"Missing ETC files")
            sys.exit(1)

        # Load datasets
        ds_mcs, ds_ar, ds_tc, ds_etc = get_datasets(dir_mcs, files_ar, files_tc, files_etc, parallel, logger)

        # Process and write output
        ds = combine_masks(ds_mcs, ds_ar, ds_tc, ds_etc, client=client, out_zarr=out_zarr, logger=logger)
        
        # Cleanup
        ds_mcs.close()
        ds_ar.close()
        ds_tc.close()
        ds_etc.close()
        if ds is not None:
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