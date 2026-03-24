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
import sys
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import healpy as hp
from easygems import healpix as egh
import time
import warnings
import json
import intake
import zarr

# Add parent directory to path to import from src
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import utility functions
from src.env_extract_utilities import (
    convert_time,
    parse_storm_ids,
    parse_pressure_levels,
    detect_pressure_units,
    normalize_pressure_levels,
    convert_w_to_omega, 
    convert_omega_to_w,
    convert_scream_zg,
    apply_model_fixes,
    parse_etc_track_file
)

# Suppress warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def extract_healpix_variable_to_latlon(
    hp_data_at_time,
    storm_lon,
    storm_lat,
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
    
    NOTE: This function now expects data already sliced to a single time step.
    Use extract_etc_2d_variable() for batched time-slice processing.
    
    Parameters
    ----------
    hp_data_at_time : xarray.DataArray
        HEALPix variable data at a single time (no 'time' dimension)
    storm_lon : float
        Storm center longitude
    storm_lat : float
        Storm center latitude
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
    # Data already sliced to single time - no time selection needed
    data_at_time = hp_data_at_time
    
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
    variable_name='var',
    unstructured_mesh=True
):
    """
    Extract 2D variable for all ETC storm positions using batched time-slice loading.
    
    This function implements the efficient batched approach:
    1. Group storms by unique timestamps
    2. For each timestamp batch, load data once: variable_data.sel(time=t).compute()
    3. Extract all storms at that timestamp from the loaded time slice
    
    This avoids redundant loading of the same time slice for multiple storms.
    
    Parameters:
    -----------
    storm_df : pd.DataFrame
        DataFrame with storm positions (storm_id, lon, lat, base_time)
    variable_data : xarray.DataArray
        Variable to extract (lazy dask array)
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
    tuple : (output_array, time_array, storm_ids, grid_ids, storm_lats, storm_lons, x_coords, y_coords, var_attrs)
    """
    # Calculate expected grid dimensions
    ny = int(2 * radius / lat_res + 1)
    nx = int(2 * radius / lon_res + 1)
    
    print(f"Extracting {variable_name} for {len(storm_df)} storm positions")
    print(f"Output grid dimensions: ({ny}, {nx})")
    print(f"Expected memory per variable: {len(storm_df) * ny * nx * 4 / 1e9:.2f} GB (float32)")
    sys.stdout.flush()
    
    # Store variable attributes from the source data
    var_attrs = dict(variable_data.attrs) if hasattr(variable_data, 'attrs') else {}
    
    # Initialize output array [time, y, x]
    ntimes = len(storm_df)
    output_array = np.full((ntimes, ny, nx), np.nan, dtype=np.float32)
    
    # Store metadata arrays (pre-allocate for efficiency)
    time_array = [None] * ntimes
    storm_ids = [None] * ntimes
    grid_ids = [None] * ntimes
    storm_lats = [None] * ntimes
    storm_lons = [None] * ntimes
    
    # Create relative coordinates once
    nx_half = nx // 2
    ny_half = ny // 2
    x_coords = np.arange(-nx_half, nx_half + 1)
    y_coords = np.arange(-ny_half, ny_half + 1)
    
    # ================================================================
    # STEP 1: Group storms by timestamp (KEY OPTIMIZATION!)
    # ================================================================
    print("Grouping storms by timestamp...")
    sys.stdout.flush()
    
    # Create mapping: timestamp -> list of (storm_index, storm_info)
    time_to_storms = {}
    for idx, row in storm_df.iterrows():
        storm_time = pd.Timestamp(row['base_time'])
        if storm_time not in time_to_storms:
            time_to_storms[storm_time] = []
        time_to_storms[storm_time].append((idx, row))
    
    unique_times = sorted(time_to_storms.keys())
    print(f"Found {len(unique_times)} unique timestamps for {len(storm_df)} storm positions")
    print(f"Average storms per timestamp: {len(storm_df) / len(unique_times):.1f}")
    sys.stdout.flush()
    
    # ================================================================
    # STEP 2: Process each timestamp batch
    # ================================================================
    start_time = time.time()
    storms_processed = 0
    
    for time_idx, storm_time in enumerate(unique_times):
        if time_idx % max(1, len(unique_times) // 20) == 0:  # Print ~20 progress updates
            elapsed = time.time() - start_time
            rate = storms_processed / elapsed if elapsed > 0 else 0
            remaining_storms = ntimes - storms_processed
            eta = remaining_storms / rate if rate > 0 else 0
            print(f"  Time batch {time_idx + 1}/{len(unique_times)}: "
                  f"{storms_processed}/{ntimes} storms ({100*storms_processed/ntimes:.1f}%) - "
                  f"{rate:.1f} storms/s - ETA: {eta/60:.1f} min")
            sys.stdout.flush()
        
        storms_at_this_time = time_to_storms[storm_time]
        
        # ================================================================
        # LOAD TIME SLICE ONCE (this is the key optimization!)
        # ================================================================
        try:
            # Select time but DON'T compute yet - keep as lazy dask array
            # IMPORTANT: .squeeze() to remove extra dimensions (e.g., pressure)
            data_at_time_lazy = variable_data.sel(time=storm_time, method='nearest').squeeze()
            
            # For efficiency: collect all pixel indices needed for all storms at this time
            # Then compute once with all pixels selected
            all_pixel_sets = []
            storm_infos = []
            
            for storm_idx, row in storms_at_this_time:
                storm_lon = row['lon']
                storm_lat = row['lat']
                
                # For structured mesh (ERA5), snap storm center to nearest HEALPix cell
                # This ensures we use actual HEALPix grid points rather than interpolating
                if not unstructured_mesh:
                    # Find nearest HEALPix cell to storm center
                    center_pix = hp.ang2pix(nside, storm_lon, storm_lat, nest=True, lonlat=True)                    
                    # Get the actual lat/lon of this HEALPix cell center
                    storm_lon_hp, storm_lat_hp = hp.pix2ang(nside, center_pix, nest=True, lonlat=True)
                    # Use HEALPix cell center for grid calculation
                    lon_center = round(storm_lon_hp / lon_res) * lon_res
                    lat_center = round(storm_lat_hp / lat_res) * lat_res
                else:
                    # For unstructured mesh, use storm position directly
                    lon_center = round(storm_lon / lon_res) * lon_res
                    lat_center = round(storm_lat / lat_res) * lat_res
                
                # Calculate grid for this storm
                lon_grid = np.arange(lon_center - radius, lon_center + radius + lon_res, lon_res)
                lat_grid = np.arange(lat_center - radius, lat_center + radius + lat_res, lat_res)
                
                # Identify valid latitudes
                valid_lat_mask = (lat_grid >= -90.0) & (lat_grid <= 90.0)
                
                if np.all(valid_lat_mask):
                    # All latitudes valid
                    pix = hp.ang2pix(nside, *np.meshgrid(lon_grid, lat_grid), 
                                    nest=True, lonlat=True)
                    storm_infos.append({
                        'storm_idx': storm_idx,
                        'row': row,
                        'pix': pix,
                        'lon_grid': lon_grid,
                        'lat_grid': lat_grid,
                        'valid_lat_mask': valid_lat_mask,
                        'needs_padding': False
                    })
                else:
                    # Some latitudes invalid - need padding
                    lat_grid_valid = lat_grid[valid_lat_mask]
                    pix = hp.ang2pix(nside, *np.meshgrid(lon_grid, lat_grid_valid), 
                                    nest=True, lonlat=True)
                    storm_infos.append({
                        'storm_idx': storm_idx,
                        'row': row,
                        'pix': pix,
                        'lon_grid': lon_grid,
                        'lat_grid': lat_grid,
                        'lat_grid_valid': lat_grid_valid,
                        'valid_lat_mask': valid_lat_mask,
                        'needs_padding': True
                    })
                
                # Collect unique pixels
                all_pixel_sets.append(set(pix.flatten()))
            
            # Get union of all pixels needed at this time
            all_pixels_needed = set()
            for pixel_set in all_pixel_sets:
                all_pixels_needed.update(pixel_set)
            all_pixels_needed = list(all_pixels_needed)
            
            # COMPUTE ONCE with all needed pixels (CRITICAL OPTIMIZATION!)
            data_at_time_subset = data_at_time_lazy.isel(cell=all_pixels_needed).compute()
            
            # Create reverse mapping: cell_id -> index in subset
            cell_to_subset_idx = {cell_id: i for i, cell_id in enumerate(all_pixels_needed)}
            
        except Exception as e:
            print(f"  ERROR: Failed to load time slice at {storm_time}: {e}")
            storms_processed += len(storms_at_this_time)
            continue
        
        # ================================================================
        # EXTRACT all storms at this timestamp from the loaded time slice
        # ================================================================
        for storm_info in storm_infos:
            storm_idx = storm_info['storm_idx']
            row = storm_info['row']
            pix = storm_info['pix']
            
            try:
                if not storm_info['needs_padding']:
                    # No padding needed - simple case
                    # Map pixel indices to subset indices
                    subset_indices = [[cell_to_subset_idx[cell] for cell in row] 
                                     for row in pix]
                    
                    # Extract from loaded subset using numpy indexing (FAST!)
                    data_gridded = data_at_time_subset.values[subset_indices]
                    
                else:
                    # Padding needed for invalid latitudes
                    lat_grid = storm_info['lat_grid']
                    lon_grid = storm_info['lon_grid']
                    valid_lat_mask = storm_info['valid_lat_mask']
                    
                    # Extract valid data
                    subset_indices = [[cell_to_subset_idx[cell] for cell in row] 
                                     for row in pix]
                    data_valid = data_at_time_subset.values[subset_indices]
                    
                    # Create full array with NaN padding
                    ny_full = len(lat_grid)
                    nx_full = len(lon_grid)
                    data_gridded = np.full((ny_full, nx_full), np.nan, dtype=data_valid.dtype)
                    
                    # Fill valid rows
                    valid_lat_indices = np.where(valid_lat_mask)[0]
                    for i, lat_idx in enumerate(valid_lat_indices):
                        data_gridded[lat_idx, :] = data_valid[i, :]
                
                # Store in output array
                output_array[storm_idx, :, :] = data_gridded
                
                # Store metadata
                time_array[storm_idx] = storm_time
                storm_ids[storm_idx] = row['storm_id']
                if unstructured_mesh:
                    grid_ids[storm_idx] = row['grid_id']
                else:
                    # For structured mesh, store as tuple (lon_id, lat_id)
                    grid_ids[storm_idx] = (row['lon_id'], row['lat_id'])
                storm_lats[storm_idx] = row['lat']
                storm_lons[storm_idx] = row['lon']
                
            except Exception as e:
                print(f"  Warning: Failed to extract storm {row['storm_id']} at {storm_time}: {e}")

                # Metadata already initialized to None, will be filled with defaults
                time_array[storm_idx] = storm_time
                storm_ids[storm_idx] = row['storm_id']
                if unstructured_mesh:
                    grid_ids[storm_idx] = row['grid_id']
                else:
                    grid_ids[storm_idx] = (row['lon_id'], row['lat_id'])
                storm_lats[storm_idx] = row['lat']
                storm_lons[storm_idx] = row['lon']
        
        storms_processed += len(storms_at_this_time)
    
    elapsed = time.time() - start_time
    print(f"Extraction complete in {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    print(f"Processing rate: {ntimes/elapsed:.1f} storms/second")
    print(f"Time slices loaded: {len(unique_times)} (vs {ntimes} in old approach)")
    print(f"Speedup from batching: {ntimes/len(unique_times):.1f}x fewer data loads")
    sys.stdout.flush()
    
    return (output_array, time_array, storm_ids, grid_ids, storm_lats, storm_lons, 
            x_coords, y_coords, var_attrs)


def save_to_zarr(output_array, time_array, storm_ids, grid_ids, storm_lats, storm_lons,
                 x_coords, y_coords, variable_name, output_path, 
                 radius, lon_res, lat_res, chunk_size=1000, unstructured_mesh=True, var_attrs=None):
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
    var_attrs : dict, optional
        Variable attributes from source data
    
    Returns:
    --------
    str : Path to saved zarr file
    """
    print(f"Saving {variable_name} to zarr format...")
    sys.stdout.flush()
    
    # Prepare variable attributes (use source attributes if available)
    if var_attrs is None:
        var_attrs = {}
    # Add extraction-specific metadata to attributes
    extraction_attrs = dict(var_attrs)  # Copy original attributes
    extraction_attrs['extraction_info'] = 'Extracted around ETC storm center'
    if 'long_name' not in extraction_attrs:
        extraction_attrs['long_name'] = f'{variable_name}'
    
    # Create xarray Dataset
    data_vars = {
        variable_name: (['time', 'y', 'x'], output_array, extraction_attrs),
        'storm_id': (['time'], np.array(storm_ids), {
            'long_name': 'Storm ID',
            'description': 'ETC storm identifier'
        }),
        'storm_lat': (['time'], np.array(storm_lats), {
            'long_name': 'Storm Center Latitude',
            'units': 'degrees_north'
        }),
        'storm_lon': (['time'], np.array(storm_lons), {
            'long_name': 'Storm Center Longitude',
            'units': 'degrees_east'
        })
    }
    
    # Add grid identifiers based on mesh type
    if unstructured_mesh:
        data_vars['grid_id'] = (['time'], np.array(grid_ids), {
            'long_name': 'HEALPix Grid ID',
            'description': 'HEALPix cell index at storm center'
        })
    else:
        # For structured mesh, extract lon_id and lat_id from tuples
        lon_ids = [gid[0] for gid in grid_ids]
        lat_ids = [gid[1] for gid in grid_ids]
        data_vars['lon_id'] = (['time'], np.array(lon_ids), {
            'long_name': 'Longitude Grid Index',
            'description': 'Longitude grid index at storm center'
        })
        data_vars['lat_id'] = (['time'], np.array(lat_ids), {
            'long_name': 'Latitude Grid Index',
            'description': 'Latitude grid index at storm center'
        })
    
    ds = xr.Dataset(
        data_vars,
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
    parser.add_argument('--current_location', default=None, 
                        help='Current location in catalog (default: None, opens catalog directly)')
    parser.add_argument('--catalog_model', required=True,  
                        help='Model name in the catalog')
    parser.add_argument('--catalog_params', default='{"zoom": 8}', 
                        help='JSON string of catalog parameters')
    parser.add_argument('--cof_mask', action='store_true', default=False,
                       help='Extract COF mask (default: False)')
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
    
    # 3D variable options
    parser.add_argument('--pressure_levels', default=None, 
                        help='Comma-separated list of pressure levels in hPa (e.g., "850,500,300")')
    
    # Vertical velocity conversion
    parser.add_argument('--convert_wa_to_omega', action='store_true',
                        help='Convert vertical velocity (wa) to pressure velocity (omega)')
    parser.add_argument('--convert_omega_to_wa', action='store_true',
                        help='Convert pressure velocity (omega) to vertical velocity (wa)')
    
    # Pre-computed variable options
    parser.add_argument('--precomputed_dir', default=None,
                        help='Directory containing pre-computed variables')
    
    # Track file format options
    parser.add_argument('--unstructured_mesh', action='store_true', default=True,
                        help='Track file is on unstructured mesh (default: True for most models). Set to False for ERA5 lat/lon grid.')
    parser.add_argument('--structured_mesh', dest='unstructured_mesh', action='store_false',
                        help='Track file is on structured lat/lon grid (ERA5)')
    
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
    
    # Parse pressure levels if provided
    pressure_levels = None
    if args.pressure_levels:
        pressure_levels = parse_pressure_levels(args.pressure_levels)
        print(f"Using pressure levels: {pressure_levels} hPa")
        sys.stdout.flush()
    
    # Open catalog and get dataset
    if not args.cof_mask:
        print(f"Opening catalog from {args.catalog_url}")
        sys.stdout.flush()
        
        # Special case for ERA5 precipitation - use IMERG instead
        if 'era5' in args.catalog_model.lower() and 'pr' in args.variables:
            print(f"Special case detected: ERA5 model with 'pr' variable")
            print(f"Loading IMERG V7 precipitation (NOT from catalog)")
            
            # Load IMERG data for precipitation
            dir_healpix = "/pscratch/sd/w/wcmca1/GPM/healpix/"
            in_basename = "IMERG_V7_"
            time_res = "6H"
            zoom = catalog_params.get('zoom', 8)
            in_zarr = f"{dir_healpix}{in_basename}{time_res}_zoom{zoom}_20190101_20211231.zarr"
            
            print(f"  Reading: {in_zarr}")
            sys.stdout.flush()
            ds = xr.open_zarr(in_zarr, consolidated=True).pipe(
                egh.attach_coords, signed_lon=True
            )
            ds = ds.assign_coords(time=convert_time(ds.time.values))
            # Rename precipitation variable to 'pr' for consistency
            ds = ds.rename({'precipitation': 'pr'})
            
            print(f"  IMERG dataset loaded successfully")
            print(f"  Variables available: {list(ds.data_vars)}")
            sys.stdout.flush()
        
        elif 'scream' in args.catalog_model.lower() and 'pr' in args.variables:
            print(f"Special case detected: SCREAM model with 'pr' variable")
            print(f"Loading SCREAM precipitation (NOT from catalog)")
            
            dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/scream/"
            in_basename = "scream_pr"
            time_res = "6h"
            zoom = catalog_params.get('zoom', 8)
            in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
            
            print(f"  Reading: {in_zarr}")
            sys.stdout.flush()
            ds = xr.open_zarr(in_zarr, consolidated=True).pipe(
                egh.attach_coords, signed_lon=True
            )
            ds = ds.assign_coords(time=convert_time(ds.time.values))
            
            print(f"  SCREAM dataset loaded successfully")
            print(f"  Variables available: {list(ds.data_vars)}")
            sys.stdout.flush()
        
        elif 'nicam_gl11' in args.catalog_model.lower() and 'pr' in args.variables:
            print(f"Special case detected: NICAM model with 'pr' variable")
            print(f"Loading NICAM precipitation (NOT from catalog)")
            
            dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/nicam_gl11/shifted/"
            in_basename = "NICAM_pr"
            time_res = "6h"
            zoom = catalog_params.get('zoom', 8)
            in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
            
            print(f"  Reading: {in_zarr}")
            sys.stdout.flush()
            ds = xr.open_zarr(in_zarr, consolidated=True).pipe(
                egh.attach_coords, signed_lon=True
            )
            ds = ds.assign_coords(time=convert_time(ds.time.values))
            
            print(f"  NICAM dataset loaded successfully")
            print(f"  Variables available: {list(ds.data_vars)}")
            sys.stdout.flush()
        
        elif 'um_glm_n2560_ral3p3' in args.catalog_model.lower() and 'pr' in args.variables:
            print(f"Special case detected: UM GLM model with 'pr' variable")
            print(f"Loading UM GLM precipitation (NOT from catalog)")
            
            dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/um_glm_n2560_RAL3p3/"
            in_basename = "um_glm_n2560_RAL3p3_pr"
            time_res = "6h"
            zoom = catalog_params.get('zoom', 8)
            in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
            
            print(f"  Reading: {in_zarr}")
            sys.stdout.flush()
            ds = xr.open_zarr(in_zarr, consolidated=True).pipe(
                egh.attach_coords, signed_lon=True
            )
            ds = ds.assign_coords(time=convert_time(ds.time.values))
            
            print(f"  UM GLM dataset loaded successfully")
            print(f"  Variables available: {list(ds.data_vars)}")
            sys.stdout.flush()
        
        elif 'casesm2_10km_nocumulus' in args.catalog_model.lower() and 'pr' in args.variables:
            print(f"Special case detected: CASESM2 model with 'pr' variable")
            print(f"Loading CASESM2 precipitation (NOT from catalog)")
            
            dir_healpix = "/pscratch/sd/w/wcmca1/hackathon/healpix/casesm2_10km_nocumulus/"
            in_basename = "casesm2_10km_nocumulus_pr"
            time_res = "6h"
            zoom = catalog_params.get('zoom', 8)
            in_zarr = f"{dir_healpix}{in_basename}{time_res}_z{zoom}.zarr"
            
            print(f"  Reading: {in_zarr}")
            sys.stdout.flush()
            ds = xr.open_zarr(in_zarr, consolidated=True).pipe(
                egh.attach_coords, signed_lon=True
            )
            ds = ds.assign_coords(time=convert_time(ds.time.values))
            
            print(f"  CASESM2 dataset loaded successfully")
            print(f"  Variables available: {list(ds.data_vars)}")
            sys.stdout.flush()
            
        else:
            # Standard catalog loading for all other cases
            if args.current_location is not None:
                # Online catalog with current location
                cat = intake.open_catalog(args.catalog_url)[args.current_location]
            else:
                # Local catalog
                cat = intake.open_catalog(args.catalog_url)
            
            print(f"Loading dataset {args.catalog_model}...")
            sys.stdout.flush()
            ds = cat[args.catalog_model](**catalog_params).to_dask().pipe(
                egh.attach_coords, signed_lon=True
            )
            ds = ds.assign_coords(time=convert_time(ds.time.values))

            # Apply model-specific fixes
            ds = apply_model_fixes(ds, args.catalog_model)
    else:
        # Read Co-occurrence Feature Mask
        cof_root_dir = "/pscratch/sd/w/wcmca1/hackathon/cof_masks/"
        if "scream" in args.catalog_model:
            source = "scream"
        elif "era5" in args.catalog_model.lower():
            source = "IMERGv7"
        else:
            source = args.catalog_model
        cof_dir = f"{cof_root_dir}{source}_cofmasks_hp8_v1.zarr"
        ds = xr.open_zarr(cof_dir, consolidated=True).pipe(
            egh.attach_coords, signed_lon=True
        )
        ds = ds.assign_coords(time=convert_time(ds.time.values))
    # import pdb; pdb.set_trace()
    
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
    print(f"Track file format: {'unstructured mesh' if args.unstructured_mesh else 'structured lat/lon grid'}")
    sys.stdout.flush()

    storm_df = parse_etc_track_file(args.trackfile, unstructured_mesh=args.unstructured_mesh)

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
        
        # Track pressure level info for output filename and variable name
        pressure_suffix = ''
        output_variable_name = variable_name  # Default: use original variable name
        
        # =================================================================
        # HANDLE VERTICAL VELOCITY CONVERSIONS
        # =================================================================
        
        # Option 1: Convert wa to omega
        if variable_name == 'wa' and args.convert_wa_to_omega:
            print("Converting vertical velocity (wa) to pressure velocity (omega)...")
            sys.stdout.flush()
            
            if pressure_levels is None:
                print("ERROR: --pressure_levels required for wa to omega conversion")
                continue
            
            # Check if required variables exist
            if 'wa' not in ds or 'ta' not in ds:
                print(f"ERROR: Variables 'wa' and 'ta' required for omega conversion")
                continue
            
            # Convert wa to omega (returns only omega variable, not full dataset)
            variable_data = convert_w_to_omega(ds, pressure_levels)
            output_variable_name = 'omega'  # Change variable name for output
            print(f"Variable data ready (omega)")
            
            # Set pressure suffix for output filename
            if len(pressure_levels) == 1:
                pressure_suffix = f"_{int(pressure_levels[0])}hPa"
                # Update output variable name
                output_variable_name = f"omega{int(pressure_levels[0])}"
            else:
                levels_str = '-'.join([str(int(p)) for p in pressure_levels])
                pressure_suffix = f"_avg{levels_str}hPa"
                output_variable_name = f"omega{levels_str}"
        
        # Option 2: Convert omega to wa
        elif variable_name == 'omega' and args.convert_omega_to_wa:
            print("Converting pressure velocity (omega) to vertical velocity (wa)...")
            sys.stdout.flush()
            
            if pressure_levels is None:
                print("ERROR: --pressure_levels required for omega to wa conversion")
                continue
            
            # Check if required variables exist
            if 'omega' not in ds or 'ta' not in ds:
                print(f"ERROR: Variables 'omega' and 'ta' required for wa conversion")
                continue
            
            # Convert omega to wa (returns only wa variable, not full dataset)
            variable_data = convert_omega_to_w(ds, pressure_levels)
            output_variable_name = 'wa'  # Change variable name for output
            print(f"Variable data ready (wa)")
            
            # Set pressure suffix for output filename
            if len(pressure_levels) == 1:
                pressure_suffix = f"_{int(pressure_levels[0])}hPa"
                # Update output variable name
                output_variable_name = f"wa{int(pressure_levels[0])}"
            else:
                levels_str = '-'.join([str(int(p)) for p in pressure_levels])
                pressure_suffix = f"_avg{levels_str}hPa"
                output_variable_name = f"wa{levels_str}"
        
        # Option 3: Standard variable from catalog
        else:
            # Check if variable exists
            if variable_name not in ds:
                print(f"ERROR: Variable '{variable_name}' not found in dataset")
                print(f"Available variables: {list(ds.data_vars)}")
                continue
            
            # Get variable data
            print(f"Loading variable data: {variable_name}")
            sys.stdout.flush()
            variable_data = ds[variable_name]

            # Handle SCREAM specific variable fixes
            if 'scream' in args.catalog_model.lower() and 'zg' in variable_name:
                # Check if ELEV exists in dataset
                if 'ELEV' in ds:
                    print("Applying SCREAM zg variable fix (add ELEV to zg)")
                    variable_data = convert_scream_zg(variable_data, ds['ELEV'])
                    print("SCREAM zg variable fix applied")
                else:
                    print("WARNING: ELEV not found in dataset - skipping SCREAM zg correction")
                sys.stdout.flush()
            
            # Check if 3D variable (has pressure dimension)
            # Only process if not already handled by conversion (conversions already handle pressure)
            if 'pressure' in variable_data.dims:
                if pressure_levels is None:
                    print(f"WARNING: Variable '{variable_name}' has pressure dimension but no pressure levels specified")
                    print(f"Use --pressure_levels to specify pressure levels (e.g., --pressure_levels 850,500,300)")
                    print(f"Skipping {variable_name}")
                    continue
                
                print(f"3D variable detected with pressure dimension")
                print(f"Requested pressure levels: {pressure_levels} hPa")
                
                # Normalize pressure levels to match dataset units
                pressure_levels_dataset, pressure_units = normalize_pressure_levels(
                    pressure_levels, variable_data.pressure
                )

                # Select specified pressure levels (using dataset units)
                variable_data = variable_data.sel(
                    pressure=pressure_levels_dataset, method='nearest'
                )
                
                # Average if multiple levels, otherwise just select single level
                if len(pressure_levels) == 1:
                    # Single level - just squeeze out the pressure dimension
                    variable_data = variable_data.squeeze('pressure', drop=True)
                    pressure_suffix = f"_{int(pressure_levels[0])}hPa"
                    # Append pressure level to variable name (e.g., hus -> hus850)
                    output_variable_name = f"{variable_name}{int(pressure_levels[0])}"
                    print(f"Selected single pressure level, new shape: {variable_data.shape}")
                    print(f"Output variable name: {output_variable_name}")
                else:
                    # Multiple levels - average them
                    variable_data = variable_data.mean(dim='pressure', keep_attrs=True)
                    levels_str = '-'.join([str(int(p)) for p in pressure_levels])
                    pressure_suffix = f"_avg{levels_str}hPa"
                    # Append averaged pressure levels to variable name (e.g., hus -> hus850-500)
                    output_variable_name = f"{variable_name}{levels_str}"
                    print(f"Averaged {len(pressure_levels)} pressure levels, new shape: {variable_data.shape}")
                    print(f"Output variable name: {output_variable_name}")

        print(f"Variable shape: {variable_data.shape}")
        print(f"Variable dimensions: {variable_data.dims}")
        sys.stdout.flush()
        
        # Extract 2D data for all storm positions using batched approach
        (output_array, time_array, storm_ids, grid_ids, storm_lats, storm_lons,
         x_coords, y_coords, var_attrs) = extract_etc_2d_variable(
            storm_df=storm_df,
            variable_data=variable_data,
            hp_grid=hp_grid,
            nside=nside,
            radius=args.radius,
            lon_res=args.lon_res,
            lat_res=args.lat_res,
            progress_freq=args.progress_freq,
            variable_name=variable_name,
            unstructured_mesh=args.unstructured_mesh
        )
        
        # Save to zarr
        start_str = args.start_date.replace('-', '') if args.start_date else "all"
        end_str = args.end_date.replace('-', '') if args.end_date else "all"
        # If filter by storm IDs is provided, take the first ID as output filename
        if storm_ids_filter:
            output_path = os.path.join(
                args.output_dir,
                f"etc_2d_{variable_name}{pressure_suffix}_track{storm_ids_filter[0]}"
            )
        else:
            output_path = os.path.join(
                args.output_dir,
                f"etc_2d_{variable_name}{pressure_suffix}_{start_str}_{end_str}"
            )
        
        zarr_path = save_to_zarr(
            output_array, time_array, storm_ids, grid_ids, storm_lats, storm_lons,
            x_coords, y_coords, output_variable_name, output_path,
            args.radius, args.lon_res, args.lat_res, args.chunk_size,
            unstructured_mesh=args.unstructured_mesh,
            var_attrs=var_attrs
        )
        
        # Print summary for this variable
        var_elapsed_time = time.time() - var_start_time
        print("\n" + "="*60)
        print(f"VARIABLE {variable_name} COMPLETED")
        if output_variable_name != variable_name:
            print(f"Saved as: {output_variable_name}")
        print("="*60)
        print(f"Output file: {zarr_path}")
        print(f"Storm points processed: {len(storm_df)}")
        print(f"Unique storms: {storm_df['storm_id'].nunique()}")
        valid_times = [t for t in time_array if t is not None]
        if valid_times:
            print(f"Time range: {min(valid_times)} to {max(valid_times)}")
        else:
            print("Time range: No valid times found")
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
