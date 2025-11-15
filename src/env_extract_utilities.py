"""
Utility functions for environmental variable extraction around ETC/AR tracks

Common functions shared between different extraction scripts:
- Time conversion
- Track file parsing
- Pressure level handling
- Model-specific fixes

Author: Zhe Feng
Last updated: November 2025
"""

import numpy as np
import pandas as pd
import xarray as xr
import sys
from datetime import datetime


def convert_time(time_array):
    """Convert cftime to standard datetime64"""
    if hasattr(time_array[0], 'year'):
        return np.array([np.datetime64(datetime(t.year, t.month, t.day, t.hour)) 
                         for t in time_array])
    return time_array


def parse_storm_ids(storm_ids_str):
    """Parse storm IDs string from command line to list of integers"""
    if storm_ids_str is None:
        return None
    try:
        return [int(sid.strip()) for sid in storm_ids_str.split(',')]
    except:
        print(f"ERROR: Could not parse storm IDs: {storm_ids_str}")
        return None


def parse_pressure_levels(pressure_str):
    """Parse pressure levels string from bash to list"""
    try:
        return [float(p.strip()) for p in pressure_str.split(',')]
    except:
        return [850, 500, 300]  # Default pressure levels


def detect_pressure_units(pressure_coord):
    """
    Detect whether pressure coordinates are in Pascals or hectopascals.
    
    Parameters:
    -----------
    pressure_coord : xarray.DataArray
        Pressure coordinate from dataset
    
    Returns:
    --------
    str : 'Pa' or 'hPa'
    
    Notes:
    ------
    Heuristic: If the median pressure value is > 2000, assume Pascals (typical range: 100000-10000 Pa)
               If the median pressure value is < 2000, assume hectopascals (typical range: 1000-100 hPa)
    """
    median_pressure = float(np.median(pressure_coord.values))
    
    if median_pressure > 2000:
        units = 'Pa'
        print(f"  Detected pressure units: Pascals (median value: {median_pressure:.1f} Pa)")
    else:
        units = 'hPa'
        print(f"  Detected pressure units: hectopascals (median value: {median_pressure:.1f} hPa)")
    
    sys.stdout.flush()
    return units


def normalize_pressure_levels(pressure_levels_hPa, dataset_pressure_coord):
    """
    Convert user-specified pressure levels (always in hPa) to match dataset units.
    
    Parameters:
    -----------
    pressure_levels_hPa : list
        Pressure levels specified by user in hPa (e.g., [850, 500, 300])
    dataset_pressure_coord : xarray.DataArray
        Pressure coordinate from the dataset
    
    Returns:
    --------
    tuple : (pressure_levels_dataset, units)
        pressure_levels_dataset : list - Pressure levels in dataset units
        units : str - Units detected ('Pa' or 'hPa')
    """
    units = detect_pressure_units(dataset_pressure_coord)
    
    if units == 'Pa':
        pressure_levels_dataset = [p * 100 for p in pressure_levels_hPa]
        print(f"  Converting pressure levels from hPa to Pa: {pressure_levels_hPa} hPa → {pressure_levels_dataset} Pa")
    else:
        pressure_levels_dataset = pressure_levels_hPa
        print(f"  Pressure levels: {pressure_levels_hPa} hPa (no conversion needed)")
    
    sys.stdout.flush()
    return pressure_levels_dataset, units


def apply_model_fixes(ds, model_name):
    """
    Apply model-specific fixes for dimension and variable names.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Input dataset
    model_name : str
        Model name (e.g., 'ifs', 'scream', 'nicam', 'icon', 'um')
    
    Returns:
    --------
    xarray.Dataset
        Fixed dataset with standardized names
    """
    print(f"\nApplying model-specific fixes for: {model_name}")
    sys.stdout.flush()
    
    # ===== FIX FOR IFS MODEL =====
    # IFS uses 'value' and 'cell' dimensions instead of standard names
    if 'value' in ds.dims and 'cell' in ds.dims:
        print("  Detected IFS model format (value/cell dimensions)")
        sys.stdout.flush()
        
        # Rename dimensions
        ds = ds.rename({'value': 'time', 'cell': 'ncells'})
        print("  Renamed: 'value' → 'time', 'cell' → 'ncells'")
        
        # IFS may have different variable names
        var_name_mapping = {
            't': 'ta',      # temperature
            'tcwv': 'prw',      # total column water vapor
            'w': 'omega',      # vertical velocity
            'q': 'hus',    # specific humidity
            'r': 'hur',     # relative humidity
            '2t': 'tas'  ,   # 2m temperature
            'sp': 'ps',     # surface pressure
            '10u': 'uas',   # 10m eastward wind
            '10v': 'vas',   # 10m northward wind
            '2d': 'tdas'    # 2m dew point temperature
        }
        
        vars_to_rename = {}
        for old_name, new_name in var_name_mapping.items():
            if old_name in ds.data_vars or old_name in ds.coords:
                vars_to_rename[old_name] = new_name
        
        if vars_to_rename:
            ds = ds.rename(vars_to_rename)
            print(f"  Renamed variables: {vars_to_rename}")
        
        sys.stdout.flush()
    
    # ===== FIX FOR PRESSURE DIMENSION NAMES =====
    # Standardize pressure dimension to 'pressure'
    if 'lev' in ds.dims:
        ds = ds.rename({'lev': 'pressure'})
        print("  Renamed dimension: 'lev' → 'pressure'")
        sys.stdout.flush()
    
    # Some models use 'level' instead of 'pressure'
    if 'level' in ds.dims:
        ds = ds.rename({'level': 'pressure'})
        print("  Renamed dimension: 'level' → 'pressure'")
        ds = ds.assign_coords(pressure=('pressure', ds.lev.values))
        ds = ds.drop_vars('lev')
        
        sys.stdout.flush()
    
    return ds


def parse_etc_track_file(file_path, unstructured_mesh=True):
    """
    Parse ETC track data from text file.
    
    Parameters:
    -----------
    file_path : str
        Path to ETC track file
    unstructured_mesh : bool
        Whether using HEALPix (True) or regular grid (False)
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: storm_id, grid_id, lon, lat, year, month, day, hour, base_time
    """
    print(f"Parsing ETC track file: {file_path}")
    sys.stdout.flush()
    
    storm_data = []
    with open(file_path, 'r') as f:
        storm_id = 0
        for line in f:
            line = line.strip()
            if line.startswith("start"):
                # New storm
                storm_id += 1
                num_timesteps, year, month, day, hour = map(int, line.split()[1:])
            else:
                # Storm details
                cols = line.split()
                if unstructured_mesh:
                    storm_data.append({
                        "storm_id": storm_id,
                        "grid_id": int(cols[0]),
                        "lon": float(cols[1]),
                        "lat": float(cols[2]),
                        "year": int(cols[-4]),
                        "month": int(cols[-3]),
                        "day": int(cols[-2]),
                        "hour": int(cols[-1]),
                        "base_time": np.datetime64(f"{int(cols[-4]):04d}-{int(cols[-3]):02d}-{int(cols[-2]):02d}T{int(cols[-1]):02d}:00")
                    })
                else:
                    storm_data.append({
                        "storm_id": storm_id,
                        "lon_id": int(cols[0]),
                        "lat_id": int(cols[1]),
                        "lon": float(cols[2]),
                        "lat": float(cols[3]),
                        "year": int(cols[-4]),
                        "month": int(cols[-3]),
                        "day": int(cols[-2]),
                        "hour": int(cols[-1]),
                        "base_time": np.datetime64(f"{int(cols[-4]):04d}-{int(cols[-3]):02d}-{int(cols[-2]):02d}T{int(cols[-1]):02d}:00")
                    })
    
    df = pd.DataFrame(storm_data)
    print(f"Parsed {len(df)} storm points from {df['storm_id'].nunique()} unique storms")
    sys.stdout.flush()
    
    return df
