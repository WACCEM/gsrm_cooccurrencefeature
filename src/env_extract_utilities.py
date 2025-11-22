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


def convert_w_to_omega(ds, pressure_levels_hPa):
    """
    Convert vertical velocity (w) to pressure velocity (omega) using ω = -ρgw
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Dataset containing 'wa' (vertical velocity) and 'ta' (temperature)
    pressure_levels_hPa : list
        List of pressure levels in hPa (will be converted to dataset units automatically)
    
    Returns:
    --------
    xarray.DataArray
        Omega variable (pressure velocity in Pa/s)
    """
    print("Converting vertical velocity (wa) to pressure velocity (omega)...")
    sys.stdout.flush()
    
    # Physical constants
    g = 9.81  # gravitational acceleration (m/s²)
    R = 287.04  # specific gas constant for dry air (J/(kg·K))
    
    # Get variables
    w = ds['wa']  # vertical velocity (m/s)
    T = ds['ta']  # temperature (K)
    
    # Normalize pressure levels to dataset units
    pressure_levels_dataset, pressure_units = normalize_pressure_levels(
        pressure_levels_hPa, w.pressure
    )

    # Create omega variable for each pressure level
    omega_levels = []
    
    for i, pressure_hPa in enumerate(pressure_levels_hPa):
        pressure_dataset = pressure_levels_dataset[i]
        
        print(f"  Processing level {pressure_hPa} hPa ({pressure_dataset} {pressure_units})...")
        sys.stdout.flush()
        
        # Select nearest pressure level
        w_level = w.sel(pressure=pressure_dataset, method='nearest')
        T_level = T.sel(pressure=pressure_dataset, method='nearest')
        
        # Convert pressure to Pa for density calculation if needed
        if pressure_units == 'hPa':
            pressure_Pa = pressure_dataset * 100
        else:
            pressure_Pa = pressure_dataset
        
        # Calculate air density: ρ = p / (R * T)
        rho = pressure_Pa / (R * T_level)
        
        # Calculate omega: ω = -ρgw
        omega_level = -rho * g * w_level
        omega_level = omega_level.assign_coords(pressure=pressure_hPa)
        
        omega_levels.append(omega_level)
    
    # Concatenate all pressure levels
    omega_combined = xr.concat(omega_levels, dim='pressure')
    
    # Add proper attributes
    omega_combined.attrs = {
        'long_name': 'Pressure velocity (omega)',
        'units': 'Pa/s',
        'description': 'Pressure velocity calculated from vertical velocity using ω = -ρgw',
        'formula': 'omega = -density * 9.81 * vertical_velocity',
        'pressure_levels_hPa': str(pressure_levels_hPa),
        'source_pressure_units': pressure_units
    }
    
    print(f"Omega conversion complete. Pressure levels: {pressure_levels_hPa} hPa")
    sys.stdout.flush()
    
    return omega_combined


def convert_omega_to_w(ds, pressure_levels_hPa):
    """
    Convert pressure velocity (omega) to vertical velocity (w) using w = -ω/(ρg)
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Dataset containing 'omega' (pressure velocity) and 'ta' (temperature)
    pressure_levels_hPa : list
        List of pressure levels in hPa (will be converted to dataset units automatically)
    
    Returns:
    --------
    xarray.DataArray
        Vertical velocity variable (m/s)
    """
    print("Converting pressure velocity (omega) to vertical velocity (wa)...")
    sys.stdout.flush()
    
    # Physical constants
    g = 9.81  # gravitational acceleration (m/s²)
    R = 287.04  # specific gas constant for dry air (J/(kg·K))
    
    # Get variables
    omega = ds['omega']  # pressure velocity (Pa/s)
    T = ds['ta']  # temperature (K)
    
    # Normalize pressure levels to dataset units
    pressure_levels_dataset, pressure_units = normalize_pressure_levels(
        pressure_levels_hPa, omega.pressure
    )
    
    # Create wa variable for each pressure level
    wa_levels = []
    
    for i, pressure_hPa in enumerate(pressure_levels_hPa):
        pressure_dataset = pressure_levels_dataset[i]
        
        print(f"  Processing level {pressure_hPa} hPa ({pressure_dataset} {pressure_units})...")
        sys.stdout.flush()
        
        # Select nearest pressure level
        omega_level = omega.sel(pressure=pressure_dataset, method='nearest')
        T_level = T.sel(pressure=pressure_dataset, method='nearest')
        
        # Convert pressure to Pa for density calculation if needed
        if pressure_units == 'hPa':
            pressure_Pa = pressure_dataset * 100
        else:
            pressure_Pa = pressure_dataset
        
        # Calculate air density: ρ = p / (R * T)
        rho = pressure_Pa / (R * T_level)
        
        # Calculate wa: w = -ω/(ρg)
        wa_level = -omega_level / (rho * g)
        wa_level = wa_level.assign_coords(pressure=pressure_hPa)
        
        wa_levels.append(wa_level)
    
    # Concatenate all pressure levels
    wa_combined = xr.concat(wa_levels, dim='pressure')
    
    # Add proper attributes
    wa_combined.attrs = {
        'long_name': 'Vertical velocity (wa)',
        'units': 'm/s',
        'description': 'Vertical velocity calculated from pressure velocity using w = -ω/(ρg)',
        'formula': 'wa = -omega / (density * 9.81)',
        'pressure_levels_hPa': str(pressure_levels_hPa),
        'source_pressure_units': pressure_units
    }
    
    print(f"Vertical velocity conversion complete. Pressure levels: {pressure_levels_hPa} hPa")
    sys.stdout.flush()
    
    return wa_combined


def compute_surface_wind_speed(ds):
    """
    Compute surface wind speed (sfcWind) from horizontal wind components.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Must contain 'uas' and 'vas' variables (eastward and northward surface winds)
    
    Returns:
    --------
    xarray.DataArray
        Surface wind speed in m/s
    """
    print("  Computing surface wind speed (sfcWind) from uas and vas...")
    sys.stdout.flush()
    
    # Check if both components are available
    if 'uas' not in ds or 'vas' not in ds:
        raise ValueError("Both 'uas' and 'vas' required to compute sfcWind")
    
    # Compute wind speed: sqrt(u^2 + v^2)
    sfcWind = np.sqrt(ds['uas']**2 + ds['vas']**2)
    
    sfcWind.attrs = {
        'long_name': 'Near-Surface Wind Speed',
        'units': 'm s-1',
        'standard_name': 'wind_speed',
        'description': 'Surface wind speed computed from eastward and northward components',
        'formula': 'sqrt(uas^2 + vas^2)'
    }
    
    print("  Surface wind speed computed")
    sys.stdout.flush()
    return sfcWind


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
    # Handle various naming conventions: lev, level, levels, plev, etc.
    
    pressure_dim_candidates = ['lev', 'level', 'levels', 'plev', 'plevs', 'p_levs']
    pressure_coord_candidates = ['lev', 'level', 'levels', 'plev', 'plevs', 'p_levs']
    
    # Find which dimension name is used
    pressure_dim = None
    for dim in pressure_dim_candidates:
        if dim in ds.dims:
            pressure_dim = dim
            break
    
    if pressure_dim and pressure_dim != 'pressure':
        # Find corresponding coordinate variable (may have different name than dimension)
        pressure_coord = None
        for coord in pressure_coord_candidates:
            if coord in ds.coords or coord in ds.data_vars:
                pressure_coord = coord
                break
        
        # If no explicit coordinate found, the dimension itself is the coordinate
        if pressure_coord is None:
            pressure_coord = pressure_dim
        
        # Rename dimension
        ds = ds.rename({pressure_dim: 'pressure'})
        print(f"  Renamed dimension: '{pressure_dim}' → 'pressure'")
        
        # Handle coordinate variable if different from dimension
        if pressure_coord != pressure_dim and pressure_coord in ds.coords:
            # Coordinate exists with different name - assign it to 'pressure' dimension
            ds = ds.assign_coords(pressure=('pressure', ds[pressure_coord].values))
            ds = ds.drop_vars(pressure_coord)
            print(f"  Assigned coordinate '{pressure_coord}' to 'pressure' dimension and dropped '{pressure_coord}'")
        elif pressure_coord != 'pressure' and pressure_coord in ds.coords:
            # Coordinate has same name as old dimension, just drop it if it still exists
            if pressure_coord in ds.coords and pressure_coord != 'pressure':
                ds = ds.drop_vars(pressure_coord)
                print(f"  Dropped redundant coordinate '{pressure_coord}'")
        
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
