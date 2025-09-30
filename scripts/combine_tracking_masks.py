import xarray as xr
import numpy as np
import pandas as pd
import cftime
import os, glob
import time
import logging
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from src.zarr_tools import write_zarr, setup_dask_client

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

    drop_var_list = ['ccs_mask', 'pr']
    rename_dict = {
        'AR_count_index': 'ar_mask',
        'TC_count_index': 'tc_mask',
        'ETC_count_index': 'etc_mask',
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

    # Write to Zarr
    if out_zarr:
        write_zarr(ds, out_zarr, client=client, logger=logger)
    
    return ds

#-------------------------------------------------------------------
def convert_to_matching_calendar(std_times, target_calendar):
    """
    Convert standard calendar (proleptic_gregorian) timestamps to match a target calendar.
    
    Args:
        std_times: array of numpy.datetime64, pandas.DatetimeIndex or pandas.Timestamp
            Timestamps with standard (proleptic_gregorian) calendar
        target_calendar: str
            Target calendar to convert to ('365_day', '360_day', 'noleap', etc.)
            
    Returns:
        cftime.datetime objects using the target calendar
    """

    
    # Initialize the appropriate cftime date type based on target calendar
    calendar_types = {
        '365_day': cftime.DatetimeNoLeap,
        'noleap': cftime.DatetimeNoLeap,
        '360_day': cftime.Datetime360Day,
        'all_leap': cftime.DatetimeAllLeap,
        'julian': cftime.DatetimeJulian,
        # Add other calendars as needed
    }
    
    if target_calendar in ['proleptic_gregorian', 'gregorian', 'standard']:
        # No conversion needed
        return std_times
    
    if target_calendar not in calendar_types:
        raise ValueError(f"Unsupported calendar: {target_calendar}")
        
    datetime_type = calendar_types[target_calendar]
    
    # Check if input is a single timestamp
    is_single_object = not hasattr(std_times, '__iter__') or isinstance(std_times, pd.Timestamp)
    
    # Convert to list for uniform processing
    times_list = [std_times] if is_single_object else std_times
    
    # Convert each timestamp to the target calendar
    converted_times = []
    for t in times_list:
        # Convert numpy.datetime64 to pandas.Timestamp which has the necessary attributes
        if isinstance(t, np.datetime64):
            ts = pd.Timestamp(t)
            converted_times.append(datetime_type(
                ts.year, ts.month, ts.day, 
                ts.hour, ts.minute, ts.second
            ))
        else:
            # For pandas.Timestamp or datetime objects that already have year, month attributes
            converted_times.append(datetime_type(
                t.year, t.month, t.day, 
                t.hour, t.minute, t.second
            ))
    
    # Return a single object or a list based on input type
    if is_single_object:
        return converted_times[0]
    else:
        return converted_times

def setup_logging():
    """
    Set the logging message level

    Args:
        None.

    Returns:
        None.
    """
    logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO)


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
    
    # Configuration parameters
    zoom = 8
    version = 'v1'
    parallel = True
    n_workers = 32
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

    source_name = "casesm2_10km_nocumulus"
    source_res = f"hp{zoom}_H"
    dir_mcs = f"/pscratch/sd/w/wcmca1/hackathon/mcs/{source_name}/mcstracking/casesm2_hrly_mcsmask_hp{zoom}_v1.zarr"
    dir_te = f"/pscratch/sd/b/beharrop/kmscale_hackathon/hackathon_pre/{source_name}_testpy/"
    basename_ar = f"AR_tracks_{source_name}_{source_res}."
    basename_tc = f"TC_tracks_{source_name}_{source_res}."
    basename_etc = f"ETC_tracks_{source_name}_{source_res}."

    # Output paths
    out_dir = "/pscratch/sd/w/wcmca1/hackathon/allmasks/"
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
            logger.info("Shutting down Dask client")
            client.close()
    
    # Log completion time
    end_time = time.time()
    elapsed_time = end_time - start_time
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    logger.info(f"Combine completed in {int(hours):02}:{int(minutes):02}:{int(seconds):02} (hh:mm:ss).")

if __name__ == "__main__":
    main()