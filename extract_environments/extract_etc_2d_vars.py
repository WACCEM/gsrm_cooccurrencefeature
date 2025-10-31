"""
Extract 2D environmental variables around ETC tracks

This script extracts 2D spatial fields (not averaged) from HEALPix model output
for each ETC track position. Unlike MCS extraction (get_env_vars.py), this:
- Extracts full 2D spatial grids centered at each storm position
- Does not average over circular areas (saves spatial structure)
- Does not extract pre-convective data
- Uses relative x/y coordinates (grid points from storm center)
- Saves data as zarr files for efficient I/O and memory management

Memory considerations:
- For 18704 storm points × 81×81 grid × 4 bytes (float32) ≈ 4.7 GB per variable
- Zarr chunking allows processing larger-than-memory datasets
- Multiple variables can be saved separately and combined later

Author: Zhe Feng (adapted from Laura Paccini's get_env_vars.py)
Last updated: October 2025
"""
import os
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import healpy as hp
from easygems import healpix as egh
from datetime import datetime
import time
import warnings
import json
import intake
import sys
import zarr

# Suppress warnings
warnings.filterwarnings("ignore", category=FutureWarning)


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


def extract_healpix_variable_to_latlon_safe(
    hp_data,
    storm_lon,
    storm_lat,
    storm_time,
    hp_grid,
    nside,
    radius=10.0,
    lon_res=0.25,
    lat_res=0.25,
    x_dimname="lon",
    y_dimname="lat",
    x_coordname="lon",
    y_coordname="lat",
    pad_invalid_latitudes=True
):
    """
    Extract HEALPix variable and remap to regular lat-lon grid centered at storm location.
    Pads output with NaN for grids with invalid latitudes beyond poles.
    
    Parameters
    ----------
    hp_data : xarray.DataArray
        HEALPix variable data with 'cell' and 'time' dimensions
    storm_lon : float
        Storm center longitude
    storm_lat : float
        Storm center latitude
    storm_time : np.datetime64
        Storm time
    hp_grid : xarray.Dataset
        HEALPix grid with 'lat' and 'lon' coordinates
    nside : int
        HEALPix nside parameter
    radius : float, optional
        Radius around storm center in degrees (default: 10.0)
    lon_res : float, optional
        Longitude resolution in degrees (default: 0.25)
    lat_res : float, optional
        Latitude resolution in degrees (default: 0.25)
    x_dimname : str, optional
        Name for longitude dimension (default: "lon")
    y_dimname : str, optional
        Name for latitude dimension (default: "lat")
    x_coordname : str, optional
        Name for longitude coordinate (default: "lon")
    y_coordname : str, optional
        Name for latitude coordinate (default: "lat")
    pad_invalid_latitudes : bool, optional
        If True, pad with NaN for invalid latitudes beyond ±90° (default: True)
    
    Returns
    -------
    xarray.DataArray
        Variable remapped to regular lat-lon grid, padded with NaN if needed
    """
    # Select data at storm time
    data_at_time = hp_data.sel(time=storm_time, method="nearest").compute()
    
    # Find the closest grid points to the storm center
    lon_center = round(storm_lon / lon_res) * lon_res
    lat_center = round(storm_lat / lat_res) * lat_res
    
    # Create the FULL lat-lon grid (including invalid latitudes)
    lon_grid_full = np.arange(lon_center - radius, lon_center + radius + lon_res, lon_res)
    lat_grid_full = np.arange(lat_center - radius, lat_center + radius + lat_res, lat_res)
    
    if pad_invalid_latitudes:
        # Identify valid latitude indices (within [-90, 90])
        valid_lat_mask = (lat_grid_full >= -90.0) & (lat_grid_full <= 90.0)
        
        # If all latitudes are valid, proceed normally
        if np.all(valid_lat_mask):
            lon_grid = lon_grid_full
            lat_grid = lat_grid_full
            
            # Get HEALPix pixel indices for the grid points
            pix = xr.DataArray(
                hp.ang2pix(nside, *np.meshgrid(lon_grid, lat_grid), nest=True, lonlat=True),
                coords=((y_dimname, lat_grid), (x_dimname, lon_grid)),
            )
            
            # Drop lat/lon coordinates from input data to avoid conflicts during isel
            data_dropped = data_at_time.drop_vars([y_coordname, x_coordname], errors='ignore')
            
            # Remap variable to lat/lon grid
            data_gridded = data_dropped.isel(cell=pix)
            
        else:
            # Some latitudes are invalid - need to pad with NaN
            # Extract only valid latitudes
            lat_grid_valid = lat_grid_full[valid_lat_mask]
            lon_grid = lon_grid_full
            
            # Get HEALPix pixel indices for VALID grid points only
            pix_valid = xr.DataArray(
                hp.ang2pix(nside, *np.meshgrid(lon_grid, lat_grid_valid), nest=True, lonlat=True),
                coords=((y_dimname, lat_grid_valid), (x_dimname, lon_grid)),
            )
            
            # Drop lat/lon coordinates from input data to avoid conflicts during isel
            data_dropped = data_at_time.drop_vars([y_coordname, x_coordname], errors='ignore')
            
            # Remap variable to lat/lon grid for valid points
            data_valid = data_dropped.isel(cell=pix_valid)
            
            # Create full output array with NaN padding
            ny_full = len(lat_grid_full)
            nx_full = len(lon_grid)
            
            # Initialize full array with NaN
            full_shape = (ny_full, nx_full)
            data_full = np.full(full_shape, np.nan, dtype=data_valid.dtype)
            
            # Find indices where to place valid data
            valid_lat_indices = np.where(valid_lat_mask)[0]
            
            # Fill in the valid data
            for i, lat_idx in enumerate(valid_lat_indices):
                data_full[lat_idx, :] = data_valid.values[i, :]
            
            # Create output DataArray with full grid
            data_gridded = xr.DataArray(
                data_full,
                coords={
                    y_dimname: lat_grid_full,
                    x_dimname: lon_grid
                },
                dims=[y_dimname, x_dimname],
                attrs=data_valid.attrs
            )
    else:
        # No padding - clip to valid latitudes only
        lat_min = max(lat_center - radius, -90.0)
        lat_max = min(lat_center + radius, 90.0)
        
        lon_grid = lon_grid_full
        lat_grid = np.arange(lat_min, lat_max + lat_res, lat_res)
        
        # Get HEALPix pixel indices for the grid points
        pix = xr.DataArray(
            hp.ang2pix(nside, *np.meshgrid(lon_grid, lat_grid), nest=True, lonlat=True),
            coords=((y_dimname, lat_grid), (x_dimname, lon_grid)),
        )
        
        # Drop lat/lon coordinates from input data to avoid conflicts during isel
        data_dropped = data_at_time.drop_vars([y_coordname, x_coordname], errors='ignore')
        
        # Remap variable to lat/lon grid
        data_gridded = data_dropped.isel(cell=pix)
    
    return data_gridded


def extract_etc_2d_variable(
    storm_df,
    variable_data,
    hp_grid,
    nside,
    radius=10.0,
    lon_res=0.25,
    lat_res=0.25,
    progress_freq=1000,
    variable_name='var'
):
    """
    Extract 2D variable for all ETC storm positions.
    
    Parameters:
    -----------
    storm_df : pd.DataFrame
        DataFrame with storm positions (storm_id, lon, lat, base_time)
    variable_data : xarray.DataArray
        Variable to extract
    hp_grid : xarray.Dataset
        HEALPix grid
    nside : int
        HEALPix nside parameter
    radius : float
        Extraction radius in degrees
    lon_res, lat_res : float
        Grid resolution
    progress_freq : int
        How often to print progress
    variable_name : str
        Variable name for progress messages
    
    Returns:
    --------
    tuple : (output_array, time_array, storm_ids, grid_ids, storm_lats, storm_lons, x_coords, y_coords)
    """
    # Calculate expected grid dimensions
    ny = int(2 * radius / lat_res + 1)
    nx = int(2 * radius / lon_res + 1)
    
    print(f"Extracting {variable_name} for {len(storm_df)} storm positions")
    print(f"Output grid dimensions: ({ny}, {nx})")
    print(f"Expected memory per variable: {len(storm_df) * ny * nx * 4 / 1e9:.2f} GB (float32)")
    sys.stdout.flush()
    
    # Initialize output array [time, y, x]
    ntimes = len(storm_df)
    output_array = np.full((ntimes, ny, nx), np.nan, dtype=np.float32)
    
    # Store metadata
    time_array = []
    storm_ids = []
    grid_ids = []
    storm_lats = []
    storm_lons = []
    
    # Create relative coordinates once
    nx_half = nx // 2
    ny_half = ny // 2
    x_coords = np.arange(-nx_half, nx_half + 1)
    y_coords = np.arange(-ny_half, ny_half + 1)
    
    # Extract data for each storm position
    start_time = time.time()
    
    for idx, (row_idx, row) in enumerate(storm_df.iterrows()):
        if idx % progress_freq == 0:
            elapsed = time.time() - start_time
            rate = idx / elapsed if elapsed > 0 else 0
            remaining = (ntimes - idx) / rate if rate > 0 else 0
            print(f"  Progress: {idx}/{ntimes} ({100*idx/ntimes:.1f}%) - "
                  f"{rate:.1f} storms/s - ETA: {remaining/60:.1f} min")
            sys.stdout.flush()
        
        storm_id = row['storm_id']
        storm_lon = row['lon']
        storm_lat = row['lat']
        storm_time = row['base_time']
        grid_id = row['grid_id']
        
        try:
            # Extract 2D data
            data_gridded = extract_healpix_variable_to_latlon_safe(
                hp_data=variable_data,
                storm_lon=storm_lon,
                storm_lat=storm_lat,
                storm_time=storm_time,
                hp_grid=hp_grid,
                nside=nside,
                radius=radius,
                lon_res=lon_res,
                lat_res=lat_res,
                pad_invalid_latitudes=True
            )
            
            # Store data
            output_array[idx, :, :] = data_gridded.values
            
        except Exception as e:
            print(f"  Warning: Failed to extract data for storm {storm_id} at time {storm_time}: {e}")
        
        # Store metadata
        time_array.append(storm_time)
        storm_ids.append(storm_id)
        grid_ids.append(grid_id)
        storm_lats.append(storm_lat)
        storm_lons.append(storm_lon)
    
    elapsed = time.time() - start_time
    print(f"Extraction complete in {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    print(f"Processing rate: {ntimes/elapsed:.1f} storms/second")
    sys.stdout.flush()
    
    return (output_array, time_array, storm_ids, grid_ids, storm_lats, storm_lons, 
            x_coords, y_coords)


def save_to_zarr(output_array, time_array, storm_ids, grid_ids, storm_lats, storm_lons,
                 x_coords, y_coords, variable_name, output_path, 
                 radius, lon_res, lat_res, chunk_size=1000):
    """
    Save extracted data to zarr format with proper chunking.
    
    Parameters:
    -----------
    output_array : np.ndarray
        3D array [time, y, x]
    time_array : list
        Time stamps
    storm_ids, grid_ids, storm_lats, storm_lons : list
        Storm metadata
    x_coords, y_coords : np.ndarray
        Relative coordinate arrays
    variable_name : str
        Variable name
    output_path : str
        Output file path (without extension)
    radius, lon_res, lat_res : float
        Grid parameters
    chunk_size : int
        Chunk size for time dimension
    
    Returns:
    --------
    str : Path to saved zarr file
    """
    print(f"Saving {variable_name} to zarr format...")
    sys.stdout.flush()
    
    # Create xarray Dataset
    ds = xr.Dataset(
        {
            variable_name: (['time', 'y', 'x'], output_array, {
                'long_name': f'{variable_name}',
                'description': 'Extracted around ETC storm center'
            }),
            'storm_id': (['time'], np.array(storm_ids), {
                'long_name': 'Storm ID',
                'description': 'ETC storm identifier'
            }),
            'grid_id': (['time'], np.array(grid_ids), {
                'long_name': 'HEALPix Grid ID',
                'description': 'HEALPix cell index at storm center'
            }),
            'storm_lat': (['time'], np.array(storm_lats), {
                'long_name': 'Storm Center Latitude',
                'units': 'degrees_north'
            }),
            'storm_lon': (['time'], np.array(storm_lons), {
                'long_name': 'Storm Center Longitude',
                'units': 'degrees_east'
            })
        },
        coords={
            'time': np.array(time_array),
            'y': (['y'], y_coords, {
                'long_name': 'y index relative to storm center',
                'units': 'grid points',
                'description': f'Latitude index relative to storm center, spacing={lat_res}°'
            }),
            'x': (['x'], x_coords, {
                'long_name': 'x index relative to storm center',
                'units': 'grid points',
                'description': f'Longitude index relative to storm center, spacing={lon_res}°'
            })
        },
        attrs={
            'description': 'ETC storm extracted 2D data with relative coordinates',
            'radius': radius,
            'lon_res': lon_res,
            'lat_res': lat_res,
            'grid_center_x': 0,
            'grid_center_y': 0,
            'creation_time': pd.Timestamp.now().isoformat()
        }
    )
    
    # Set up chunking - chunk along time dimension for efficient time-slicing
    # Keep spatial dimensions unchunked for fast spatial access
    encoding = {
        variable_name: {
            'chunks': (chunk_size, len(y_coords), len(x_coords)),
            'compressor': zarr.Blosc(cname='zstd', clevel=3, shuffle=2)
        }
    }
    
    zarr_path = f"{output_path}.zarr"
    
    # Remove existing zarr if it exists
    if os.path.exists(zarr_path):
        import shutil
        shutil.rmtree(zarr_path)
    
    # Write to zarr
    ds.to_zarr(zarr_path, mode='w', encoding=encoding)
    
    print(f"Saved to: {zarr_path}")
    
    # Print file size
    total_size = sum(
        os.path.getsize(os.path.join(dirpath, filename))
        for dirpath, dirnames, filenames in os.walk(zarr_path)
        for filename in filenames
    )
    print(f"Total size: {total_size / 1e9:.2f} GB")
    sys.stdout.flush()
    
    return zarr_path


def main():
    parser = argparse.ArgumentParser(
        description='Extract 2D environmental variables around ETC tracks.'
    )
    
    # Input/output options
    parser.add_argument('--catalog_url', 
                        default="https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml",
                        help='URL of the intake catalog')
    parser.add_argument('--current_location', default="NERSC", 
                        help='Current location in catalog')
    parser.add_argument('--catalog_model', required=True,  
                        help='Model name in the catalog')
    parser.add_argument('--catalog_params', default='{"zoom": 8}', 
                        help='JSON string of catalog parameters')
    parser.add_argument('--trackfile', required=True, 
                        help='Path to ETC track file (.txt)')
    parser.add_argument('--output_dir', required=True, 
                        help='Output directory for results')
    parser.add_argument('--variable', '--variables', nargs='+', required=True, 
                        dest='variables',
                        help='Variable(s) to extract from dataset (can specify multiple)')
    
    # Date filtering options
    parser.add_argument('--start_date', help='Start date for filtering (YYYY-MM-DD)')
    parser.add_argument('--end_date', help='End date for filtering (YYYY-MM-DD)')
    
    # Spatial filtering options
    parser.add_argument('--min_lon', type=float, default=None, help='Minimum longitude')
    parser.add_argument('--max_lon', type=float, default=None, help='Maximum longitude')
    parser.add_argument('--min_lat', type=float, default=None, help='Minimum latitude')
    parser.add_argument('--max_lat', type=float, default=None, help='Maximum latitude')
    
    # Storm filtering (for testing)
    parser.add_argument('--storm_ids', default=None,
                        help='Comma-separated list of storm IDs to process (for testing)')
    
    # Processing options
    parser.add_argument('--radius', type=float, default=10.0,
                        help='Extraction radius in degrees (default: 10.0)')
    parser.add_argument('--lon_res', type=float, default=0.25,
                        help='Longitude resolution in degrees (default: 0.25)')
    parser.add_argument('--lat_res', type=float, default=0.25,
                        help='Latitude resolution in degrees (default: 0.25)')
    parser.add_argument('--chunk_size', type=int, default=1000,
                        help='Chunk size for time dimension in zarr (default: 1000)')
    parser.add_argument('--progress_freq', type=int, default=1000,
                        help='How often to print progress (default: 1000)')
    
    args = parser.parse_args()
    
    # Start timing
    total_start_time = time.time()
    
    print("="*60)
    print(f"ETC 2D Variable Extraction")
    print("="*60)
    print(f"Model: {args.catalog_model}")
    print(f"Track file: {args.trackfile}")
    print(f"Output directory: {args.output_dir}")
    print(f"Variables: {', '.join(args.variables)}")
    print(f"Extraction radius: {args.radius}°")
    print(f"Grid resolution: {args.lon_res}° × {args.lat_res}°")
    print("="*60)
    sys.stdout.flush()
    
    # Parse catalog parameters
    try:
        catalog_params = json.loads(args.catalog_params)
        zoom_level = catalog_params.get('zoom', 'unknown')
    except:
        catalog_params = {'zoom': 8}
        zoom_level = 8
    
    # Parse storm IDs if provided
    storm_ids_filter = parse_storm_ids(args.storm_ids)
    if storm_ids_filter:
        print(f"Testing mode: Processing only storm IDs: {storm_ids_filter}")
        sys.stdout.flush()
    
    # Open catalog and get dataset
    print(f"Opening catalog from {args.catalog_url}")
    sys.stdout.flush()
    cat = intake.open_catalog(args.catalog_url)[args.current_location]
    
    print(f"Loading dataset {args.catalog_model}...")
    sys.stdout.flush()
    ds = cat[args.catalog_model](**catalog_params).to_dask().pipe(
        egh.attach_coords, signed_lon=True
    )
    ds = ds.assign_coords(time=convert_time(ds.time.values))
    
    # Get HEALPix grid
    print("Computing HEALPix grid...")
    sys.stdout.flush()
    hp_grid = ds[['lat', 'lon']].compute()
    nside = egh.get_nside(hp_grid)
    print(f"Grid: nside={nside}")
    sys.stdout.flush()
    
    # =================================================================
    # LOAD ETC TRACK DATA ONCE (for all variables!)
    # =================================================================
    print(f"Loading ETC track data from {args.trackfile}")
    sys.stdout.flush()
    
    storm_df = parse_etc_track_file(args.trackfile, unstructured_mesh=True)
    
    # Filter by storm IDs if provided (for testing)
    if storm_ids_filter:
        storm_df = storm_df[storm_df['storm_id'].isin(storm_ids_filter)]
        print(f"Filtered to {len(storm_df)} storm points from {storm_df['storm_id'].nunique()} storms")
    
    # Apply date filtering
    if args.start_date and args.end_date:
        start_date = pd.Timestamp(args.start_date)
        end_date = pd.Timestamp(args.end_date)
        date_filter = (
            (pd.to_datetime(storm_df['base_time']) >= start_date) & 
            (pd.to_datetime(storm_df['base_time']) <= end_date)
        )
        storm_df = storm_df[date_filter]
        print(f"After date filtering: {len(storm_df)} storm points")
    
    # Apply spatial filtering
    storm_df = storm_df.dropna(subset=['lat', 'lon'])
    print(f"After dropping NaN positions: {len(storm_df)} valid storm points")
    
    # Set spatial bounds
    if args.min_lat is None:
        min_lat = -90 + args.radius
    else:
        min_lat = args.min_lat
    
    if args.max_lat is None:
        max_lat = 90 - args.radius
    else:
        max_lat = args.max_lat
    
    spatial_filter_conditions = [
        storm_df['lat'].between(min_lat, max_lat),
        storm_df['lat'].notna(),
        storm_df['lon'].notna()
    ]
    
    if args.min_lon is not None and args.max_lon is not None:
        if args.min_lon > args.max_lon:
            spatial_filter_conditions.append(
                (storm_df['lon'] >= args.min_lon) | (storm_df['lon'] <= args.max_lon)
            )
        else:
            spatial_filter_conditions.append(
                storm_df['lon'].between(args.min_lon, args.max_lon)
            )
    
    storm_df = storm_df[np.logical_and.reduce(spatial_filter_conditions)].copy()
    storm_df = storm_df.reset_index(drop=True)
    print(f"After spatial filtering: {len(storm_df)} storm points")
    print(f"Number of unique storms: {storm_df['storm_id'].nunique()}")
    print(f"This dataset will be reused for all {len(args.variables)} variable(s)")
    sys.stdout.flush()
    
    if len(storm_df) == 0:
        print("ERROR: No storm points remain after filtering!")
        return
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # =================================================================
    # PROCESS EACH VARIABLE
    # =================================================================
    for var_idx, variable_name in enumerate(args.variables):
        var_start_time = time.time()
        
        print("\n" + "="*60)
        print(f"Processing variable {var_idx + 1}/{len(args.variables)}: {variable_name}")
        print("="*60)
        sys.stdout.flush()
        
        # Check if variable exists
        if variable_name not in ds:
            print(f"ERROR: Variable '{variable_name}' not found in dataset")
            print(f"Available variables: {list(ds.data_vars)}")
            continue
        
        # Get variable data
        print(f"Loading variable data: {variable_name}")
        sys.stdout.flush()
        variable_data = ds[variable_name]
        
        # Check if 3D variable (has pressure dimension)
        if 'pressure' in variable_data.dims:
            print(f"WARNING: Variable '{variable_name}' has pressure dimension")
            print(f"This script extracts 2D variables. Use appropriate pressure level selection.")
            print(f"Skipping {variable_name}")
            continue
        
        print(f"Variable shape: {variable_data.shape}")
        print(f"Variable dimensions: {variable_data.dims}")
        sys.stdout.flush()
        
        # Extract 2D data for all storm positions
        (output_array, time_array, storm_ids, grid_ids, storm_lats, storm_lons,
         x_coords, y_coords) = extract_etc_2d_variable(
            storm_df=storm_df,
            variable_data=variable_data,
            hp_grid=hp_grid,
            nside=nside,
            radius=args.radius,
            lon_res=args.lon_res,
            lat_res=args.lat_res,
            progress_freq=args.progress_freq,
            variable_name=variable_name
        )
        
        # Save to zarr
        start_str = args.start_date.replace('-', '') if args.start_date else "all"
        end_str = args.end_date.replace('-', '') if args.end_date else "all"
        output_path = os.path.join(
            args.output_dir,
            f"etc_2d_{variable_name}_{start_str}_{end_str}"
        )
        
        zarr_path = save_to_zarr(
            output_array, time_array, storm_ids, grid_ids, storm_lats, storm_lons,
            x_coords, y_coords, variable_name, output_path,
            args.radius, args.lon_res, args.lat_res, args.chunk_size
        )
        
        # Print summary for this variable
        var_elapsed_time = time.time() - var_start_time
        print("\n" + "="*60)
        print(f"VARIABLE {variable_name} COMPLETED")
        print("="*60)
        print(f"Output file: {zarr_path}")
        print(f"Storm points processed: {len(storm_df)}")
        print(f"Unique storms: {storm_df['storm_id'].nunique()}")
        print(f"Time range: {min(time_array)} to {max(time_array)}")
        print(f"NaN values: {np.isnan(output_array).sum()} / {output_array.size} "
              f"({100*np.isnan(output_array).sum()/output_array.size:.2f}%)")
        print(f"Total time: {var_elapsed_time:.2f} seconds ({var_elapsed_time/60:.2f} minutes)")
        print("="*60)
        sys.stdout.flush()
    
    # Print final summary
    total_elapsed_time = time.time() - total_start_time
    print("\n" + "="*60)
    print("ALL VARIABLES COMPLETED")
    print("="*60)
    print(f"Processed {len(args.variables)} variable(s)")
    print(f"Total time: {total_elapsed_time:.2f} seconds ({total_elapsed_time/60:.2f} minutes)")
    if len(args.variables) > 0:
        print(f"Average per variable: {total_elapsed_time/len(args.variables):.2f} seconds")
    print("="*60)


if __name__ == "__main__":
    main()
