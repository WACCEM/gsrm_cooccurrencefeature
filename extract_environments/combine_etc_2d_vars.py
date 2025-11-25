"""
Combine individual extracted ETC 2D variable zarr files into a single multi-variable zarr file.

This script:
1. Scans a directory for etc_2d_*.zarr files
2. Groups files by their suffix pattern (e.g., track1014, all_all)
3. Combines variables from the same group into a single zarr file
4. Applies variable unit standardization based on source model/dataset
5. Optionally adds ETC COF (Co-Occurrence Feature) overlap data from parquet files:
   - Adds overlap_flag (0=isolated, 1=MCS only, 2=AR only, 3=MCS+AR)
   - Adds AR and MCS track IDs for overlapping features
   - Adds ETC center coordinates (cof_lat, cof_lon) from tracking data
   - Enabled by default, can be disabled with --no-cof-data flag
6. Preserves all metadata and coordinates

Author: Zhe Feng
Last updated: November 2025
"""
import os
import argparse
import xarray as xr
import zarr
import numpy as np
from pathlib import Path
import re
import sys
from collections import defaultdict
import pandas as pd

# Import variable scaling configuration
try:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
    from variable_scaling_config import VARIABLE_SCALING, get_scaling_info, rename_variables
except ImportError:
    # If import fails, define minimal config inline
    print("WARNING: Could not import variable_scaling_config, using minimal inline config")
    VARIABLE_SCALING = {
        'scream': {'pr': (3600000.0, 'mm h-1', 'Convert kg m-2 s-1 to mm/h'), 'zg*': (1.0/9.81, 'm', 'Convert geopotential')},
        'era5': {'pr': (1.0, 'mm h-1', 'Already in mm/h'), 'zg*': (1.0, 'm', 'Already in geopotential height')},
        'nicam': {'pr': (3600.0, 'mm h-1', 'Convert kg m-2 s-1 to mm/h'), 'zg*': (1.0/9.81, 'm', 'Convert geopotential')},
        'icon': {'pr': (3600.0, 'mm h-1', 'Convert kg m-2 s-1 to mm/h'), 'zg*': (1.0/9.81, 'm', 'Convert geopotential')},
        'cesm2': {'pr': (3600.0, 'mm h-1', 'Convert kg m-2 s-1 to mm/h'), 'zg*': (1.0, 'm', 'Already in geopotential height')},
        'um': {'pr': (3600.0, 'mm h-1', 'Convert kg m-2 s-1 to mm/h'), 'zg*': (1.0, 'm', 'Already in geopotential height')},
    }
    
    def get_scaling_info(source, variable):
        """Fallback pattern matching implementation."""
        if source not in VARIABLE_SCALING:
            return None
        config = VARIABLE_SCALING[source]
        # Try exact match first
        if variable in config:
            return config[variable]
        # Try pattern matching
        for pattern, scaling_info in config.items():
            if pattern.endswith('*') and variable.startswith(pattern[:-1]):
                return scaling_info
        return None


def standardize_variable_units(ds, source):
    """
    Standardize variable units according to source-specific scaling configuration.
    
    This function applies scaling factors to variables to convert them to 
    standardized units across different model sources. For example:
    - Precipitation: converted to mm/h
    - Geopotential: converted to geopotential height in meters
    
    Supports wildcard patterns (e.g., "zg*" matches all variables starting with "zg").
    Exact matches take precedence over pattern matches.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Dataset containing variables to scale
    source : str
        Source identifier (e.g., 'scream', 'era5', 'nicam')
    
    Returns:
    --------
    ds : xarray.Dataset
        Dataset with scaled variables and updated attributes
    """
    if source not in VARIABLE_SCALING:
        print(f"  No scaling configuration for source '{source}' - skipping unit standardization")
        return ds
    
    scaled_vars = []
    
    print(f"\n  Standardizing variable units for source: {source}")
    
    # Iterate through all data variables in the dataset
    for var_name in ds.data_vars:
        # Use get_scaling_info which handles both exact and pattern matching
        scaling_info = get_scaling_info(source, var_name)
        
        if scaling_info:
            scale_factor, target_units, description = scaling_info
            
            # Get original units
            original_units = ds[var_name].attrs.get('units', 'unknown')
            
            # Apply scaling
            if scale_factor != 1.0:
                ds[var_name] = ds[var_name] * scale_factor
                print(f"    {var_name}: {original_units} → {target_units} (×{scale_factor:.6g})")
            else:
                print(f"    {var_name}: {target_units} (no scaling needed)")
            
            # Update attributes
            ds[var_name].attrs['units'] = target_units
            ds[var_name].attrs['scaling_applied'] = description
            ds[var_name].attrs['original_units'] = original_units
            ds[var_name].attrs['scale_factor_applied'] = scale_factor
            
            scaled_vars.append(var_name)
    
    if scaled_vars:
        print(f"  Standardized {len(scaled_vars)} variable(s): {', '.join(scaled_vars)}")
    else:
        print(f"  No variables required scaling")
    
    return ds


def parse_filename(filename):
    """
    Parse zarr filename to extract variable name and suffix.
    
    Parameters:
    -----------
    filename : str
        Filename like 'etc_2d_va850_track1014.zarr'
    
    Returns:
    --------
    tuple : (variable_name, suffix)
        e.g., ('va850', 'track1014')
    """
    # Remove .zarr extension
    basename = filename.replace('.zarr', '')
    
    # Pattern: etc_2d_{variable}_{suffix}
    # Handle various patterns like track1014, all_all, YYYYMMDD_YYYYMMDD
    match = re.match(r'etc_2d_(.+?)_(track\d+|all_all|\d{8}_\d{8})', basename)
    
    if match:
        variable_name = match.group(1)
        suffix = match.group(2)
        return variable_name, suffix
    else:
        # If no standard suffix, treat everything after etc_2d_ as variable
        variable_name = basename.replace('etc_2d_', '')
        return variable_name, 'unknown'


def group_zarr_files(directory):
    """
    Group zarr files by their suffix pattern.
    
    Parameters:
    -----------
    directory : str
        Directory containing zarr files
    
    Returns:
    --------
    dict : {suffix: [(variable_name, filepath), ...]}
    """
    groups = defaultdict(list)
    
    # Find all .zarr directories
    zarr_files = [f for f in os.listdir(directory) 
                  if f.startswith('etc_2d_') and f.endswith('.zarr') 
                  and os.path.isdir(os.path.join(directory, f))]
    
    print(f"Found {len(zarr_files)} zarr files")
    
    for filename in zarr_files:
        variable_name, suffix = parse_filename(filename)
        filepath = os.path.join(directory, filename)
        groups[suffix].append((variable_name, filepath))
    
    return groups


def add_cof_data_to_combined(combined_ds, source, etc_path=None):
    """
    Add ETC COF overlap information to the combined dataset.
    
    This function loads ETC COF parquet data and adds overlap flags and track IDs
    to the combined zarr dataset by matching storm_id and time.
    
    Parameters:
    -----------
    combined_ds : xarray.Dataset
        Combined dataset with storm_id and time coordinates
    source : str
        Source identifier (e.g., 'era5', 'scream', 'nicam_gl11')
    etc_path : str, optional
        Path to ETC COF parquet files (default: /pscratch/sd/w/wcmca1/hackathon/etc_tracks/)
    
    Returns:
    --------
    combined_ds : xarray.Dataset
        Dataset with added COF variables: overlap_flag, ar_tracks_str, mcs_tracks_str, cof_lat, cof_lon
    """
    if etc_path is None:
        etc_path = '/pscratch/sd/w/wcmca1/hackathon/etc_tracks/'
    
    etc_file = f"{etc_path}/{source}_etc_cof_data.parquet"
    
    print(f"\n  Loading ETC COF data from: {etc_file}")
    
    if not os.path.isfile(etc_file):
        print(f"  WARNING: ETC COF file not found: {etc_file}")
        print(f"  Skipping COF data integration")
        return combined_ds
    
    # Load COF data
    etc_df = pd.read_parquet(etc_file)
    print(f"    Loaded {len(etc_df)} COF records")
    print(f"    Unique storms: {etc_df['storm_id'].nunique()}")
    
    # Ensure base_time in DataFrame matches the time coordinate
    etc_df['time'] = pd.to_datetime(etc_df['base_time'])
    
    # Create arrays for the new variables matching dataset's time dimension
    n_times = len(combined_ds.time)
    
    # Initialize arrays
    overlap_flag_array = np.full(n_times, np.nan)
    ar_tracks_list = [[] for _ in range(n_times)]
    mcs_tracks_list = [[] for _ in range(n_times)]
    cof_lat_array = np.full(n_times, np.nan)
    cof_lon_array = np.full(n_times, np.nan)
    
    # Create a lookup dictionary for fast access
    # Key: (storm_id, time), Value: (overlap_flag, ar_tracks, mcs_tracks, lat, lon)
    etc_lookup = {}
    for idx, row in etc_df.iterrows():
        key = (row['storm_id'], pd.Timestamp(row['time']))
        etc_lookup[key] = (row['overlap_flag'], row['ar_tracks'], row['mcs_tracks'], 
                           row['lat'], row['lon'])
    
    print(f"    Created lookup dictionary with {len(etc_lookup)} entries")
    
    # Match dataset time points with DataFrame
    for i, (time_val, storm_id_val) in enumerate(zip(combined_ds.time.values, combined_ds.storm_id.values)):
        time_key = pd.Timestamp(time_val)
        key = (int(storm_id_val), time_key)
        
        if key in etc_lookup:
            overlap_flag_array[i] = etc_lookup[key][0]
            ar_tracks_list[i] = etc_lookup[key][1]
            mcs_tracks_list[i] = etc_lookup[key][2]
            cof_lat_array[i] = etc_lookup[key][3]
            cof_lon_array[i] = etc_lookup[key][4]
    
    # Check matching statistics
    matched = np.sum(~np.isnan(overlap_flag_array))
    print(f"    Matched {matched}/{n_times} time points ({100*matched/n_times:.1f}%)")
    
    print(f"    Overlap flag distribution:")
    for flag in [0, 1, 2, 3]:
        count = np.sum(overlap_flag_array == flag)
        if matched > 0:
            print(f"      Flag {flag}: {count} ({100*count/matched:.1f}% of matched)")
        else:
            print(f"      Flag {flag}: {count}")
    
    # Add the new variables to the dataset
    combined_ds['overlap_flag'] = (['time'], overlap_flag_array, {
        'long_name': 'Co-occurrence overlap flag',
        'description': '0=isolated, 1=MCS only, 2=AR only, 3=MCS+AR',
        'units': '1'
    })
    
    # For ar_tracks and mcs_tracks, we'll store them as string representations
    ar_tracks_str = [str(tracks) if tracks is not None and len(tracks) > 0 else '[]' 
                     for tracks in ar_tracks_list]
    mcs_tracks_str = [str(tracks) if tracks is not None and len(tracks) > 0 else '[]' 
                      for tracks in mcs_tracks_list]
    
    combined_ds['ar_tracks_str'] = (['time'], ar_tracks_str, {
        'long_name': 'AR track IDs',
        'description': 'List of AR tracks overlapping with ETC (string representation)',
        'units': '1'
    })
    
    combined_ds['mcs_tracks_str'] = (['time'], mcs_tracks_str, {
        'long_name': 'MCS track IDs',
        'description': 'List of MCS tracks overlapping with ETC (string representation)',
        'units': '1'
    })
    
    combined_ds['cof_lat'] = (['time'], cof_lat_array, {
        'long_name': 'ETC center latitude from COF data',
        'description': 'Latitude of ETC center from co-occurrence feature tracking',
        'units': 'degrees_north'
    })
    
    combined_ds['cof_lon'] = (['time'], cof_lon_array, {
        'long_name': 'ETC center longitude from COF data',
        'description': 'Longitude of ETC center from co-occurrence feature tracking',
        'units': 'degrees_east'
    })
    
    print(f"    Added COF variables: overlap_flag, ar_tracks_str, mcs_tracks_str, cof_lat, cof_lon")
    
    return combined_ds


def combine_zarr_files(file_list, output_path, suffix, source=None, chunk_size=1000, add_cof_data=True, etc_path=None):
    """
    Combine multiple zarr files into a single multi-variable zarr file.
    
    Parameters:
    -----------
    file_list : list of tuples
        [(variable_name, filepath), ...]
    output_path : str
        Output zarr file path
    suffix : str
        Suffix identifier (for logging)
    source : str, optional
        Source identifier for variable scaling (e.g., 'scream', 'era5')
    chunk_size : int
        Chunk size for time dimension
    add_cof_data : bool
        Whether to add COF overlap data to the combined file (default: True)
    etc_path : str, optional
        Path to ETC COF parquet files (default: /pscratch/sd/w/wcmca1/hackathon/etc_tracks/)
    """
    print(f"\nCombining {len(file_list)} variables for suffix '{suffix}':")
    
    # Load all datasets
    datasets = {}
    reference_ds = None
    
    for var_name, filepath in sorted(file_list):
        print(f"  Loading {var_name} from {os.path.basename(filepath)}")
        sys.stdout.flush()
        
        try:
            ds = xr.open_zarr(filepath)
            
            # Use first dataset as reference for coordinates and metadata
            if reference_ds is None:
                reference_ds = ds
            
            # Check consistency
            if not np.array_equal(ds['time'].values, reference_ds['time'].values):
                print(f"    WARNING: Time coordinates differ for {var_name}")
            
            if not np.array_equal(ds['storm_id'].values, reference_ds['storm_id'].values):
                print(f"    WARNING: Storm IDs differ for {var_name}")
            
            # Extract the data variable (should be the only one besides metadata)
            # Handle both unstructured mesh (grid_id) and structured mesh (lon_id, lat_id)
            data_vars = [v for v in ds.data_vars if v not in 
                        ['storm_id', 'grid_id', 'lon_id', 'lat_id', 'storm_lat', 'storm_lon']]
            
            if len(data_vars) == 1:
                actual_var_name = data_vars[0]
                datasets[var_name] = ds[actual_var_name]
            else:
                print(f"    WARNING: Expected 1 data variable, found {len(data_vars)}: {data_vars}")
                if len(data_vars) > 0:
                    datasets[var_name] = ds[data_vars[0]]
        
        except Exception as e:
            print(f"    ERROR loading {var_name}: {e}")
            continue
    
    if not datasets:
        print(f"ERROR: No datasets loaded successfully for suffix '{suffix}'")
        return None
    
    print(f"\nSuccessfully loaded {len(datasets)} variables")
    print(f"Creating combined dataset...")
    sys.stdout.flush()
    
    # Determine if this is unstructured mesh (grid_id) or structured mesh (lon_id, lat_id)
    has_grid_id = 'grid_id' in reference_ds
    has_lon_lat_id = 'lon_id' in reference_ds and 'lat_id' in reference_ds
    
    # Create metadata dictionary based on grid type
    metadata_vars = {
        'storm_id': reference_ds['storm_id'],
        'storm_lat': reference_ds['storm_lat'],
        'storm_lon': reference_ds['storm_lon']
    }
    
    if has_grid_id:
        metadata_vars['grid_id'] = reference_ds['grid_id']
        print(f"Grid type: Unstructured mesh (HEALPix)")
    elif has_lon_lat_id:
        metadata_vars['lon_id'] = reference_ds['lon_id']
        metadata_vars['lat_id'] = reference_ds['lat_id']
        print(f"Grid type: Structured lat/lon grid")
    else:
        print(f"WARNING: No grid identifiers found (grid_id or lon_id/lat_id)")
    
    # Create combined dataset
    combined_ds = xr.Dataset(
        data_vars={
            **{var_name: (['time', 'y', 'x'], data.values, data.attrs) 
               for var_name, data in datasets.items()},
            **metadata_vars
        },
        coords={
            'time': reference_ds['time'],
            'y': reference_ds['y'],
            'x': reference_ds['x']
        },
        attrs={
            **reference_ds.attrs,
            'combined_variables': ', '.join(sorted(datasets.keys())),
            'num_variables': len(datasets),
            'combination_date': np.datetime64('now').astype(str)
        }
    )
    
    print(f"Combined dataset shape: {combined_ds.dims}")
    print(f"Variables: {list(datasets.keys())}")
    sys.stdout.flush()
    
    # Apply variable renaming if source is specified
    if source:
        try:
            combined_ds, renamed_vars = rename_variables(combined_ds, source)
            if renamed_vars:
                print(f"  Renamed {len(renamed_vars)} variable(s)")
        except NameError:
            # rename_variables not available in fallback mode
            print("\n  Variable renaming not available (using fallback config)")
    
    # Apply variable unit standardization if source is specified
    if source:
        combined_ds = standardize_variable_units(combined_ds, source)
    else:
        print("\n  No source specified - skipping variable renaming and unit standardization")
    
    # Add COF overlap data if requested
    if add_cof_data and source:
        combined_ds = add_cof_data_to_combined(combined_ds, source, etc_path)
    elif add_cof_data and not source:
        print("\n  WARNING: Cannot add COF data without source specification - skipping")
    
    # Set up chunking - same as individual files
    ny = len(combined_ds['y'])
    nx = len(combined_ds['x'])
    
    encoding = {}
    for var_name in combined_ds.keys():
        encoding[var_name] = {
            'chunks': (chunk_size, ny, nx),
            'compressor': zarr.Blosc(cname='zstd', clevel=3, shuffle=2)
        }
    
    # Remove existing output if it exists
    if os.path.exists(output_path):
        print(f"Removing existing output: {output_path}")
        import shutil
        shutil.rmtree(output_path)

    # Write to zarr
    print(f"Writing combined dataset to: {output_path}")
    sys.stdout.flush()
    
    combined_ds.to_zarr(output_path, mode='w', encoding=encoding)
    
    # Print file size
    total_size = sum(
        os.path.getsize(os.path.join(dirpath, filename))
        for dirpath, dirnames, filenames in os.walk(output_path)
        for filename in filenames
    )
    print(f"Combined file size: {total_size / 1e9:.3f} GB")
    
    return output_path


def verify_combined_file(filepath):
    """
    Verify the combined zarr file can be opened and print summary.
    
    Parameters:
    -----------
    filepath : str
        Path to combined zarr file
    """
    print(f"\nVerifying combined file: {filepath}")
    sys.stdout.flush()
    
    try:
        ds = xr.open_zarr(filepath)
        
        print(f"  Dimensions: {dict(ds.dims)}")
        print(f"  Coordinates: {list(ds.coords)}")
        
        # Determine grid type
        if 'grid_id' in ds:
            print(f"  Grid type: Unstructured mesh (HEALPix)")
            metadata_vars = ['storm_id', 'grid_id', 'storm_lat', 'storm_lon']
        elif 'lon_id' in ds and 'lat_id' in ds:
            print(f"  Grid type: Structured lat/lon grid")
            metadata_vars = ['storm_id', 'lon_id', 'lat_id', 'storm_lat', 'storm_lon']
        else:
            print(f"  Grid type: Unknown")
            metadata_vars = ['storm_id', 'storm_lat', 'storm_lon']
        
        data_vars = [v for v in ds.data_vars if v not in metadata_vars]
        print(f"  Data variables ({len(data_vars)}): {data_vars}")
        
        print(f"  Time range: {ds['time'].values[0]} to {ds['time'].values[-1]}")
        print(f"  Number of storms: {len(ds['time'])}")
        print(f"  Unique storm IDs: {len(np.unique(ds['storm_id'].values))}")
        
        # Check for NaN values in each variable
        print(f"\n  NaN statistics:")
        for var_name in data_vars:
            var_data = ds[var_name].values
            nan_count = np.isnan(var_data).sum()
            nan_pct = 100 * nan_count / var_data.size
            print(f"    {var_name}: {nan_count} / {var_data.size} ({nan_pct:.2f}%)")
        
        print(f"\n  Verification PASSED ✓")
        return True
    
    except Exception as e:
        print(f"  ERROR during verification: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Combine individual ETC 2D variable zarr files into multi-variable zarr files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Combine all zarr files for ETC data (includes COF overlap data by default)
  python combine_etc_2d_vars.py --source era5
  
  # Combine without COF overlap data
  python combine_etc_2d_vars.py --source era5 --no-cof-data
  
  # Specify custom COF data path
  python combine_etc_2d_vars.py --source era5 --etc-path /custom/path/to/etc_tracks/
  
  # Specify optional input directory for SCREAM data
  python combine_etc_2d_vars.py --input_dir /path/to/zarr/files --source scream
  
  # Specify optional output directory and source
  python combine_etc_2d_vars.py --input_dir /path/to/zarr/files --output_dir /path/to/output --source era5
  
  # Process only specific suffix pattern with source
  python combine_etc_2d_vars.py --input_dir /path/to/zarr/files --suffix track1014 --source nicam_gl11
  
  # Custom output filename prefix
  python combine_etc_2d_vars.py --input_dir /path/to/zarr/files --output_prefix etc_combined --source icon_d3hp003
        """
    )
    
    parser.add_argument('--source', required=True,
                        choices=['scream', 'era5', 'nicam_gl11', 'icon_d3hp003', 'casesm2_10km_nocumulus', 'um_glm_n2560_RAL3p3'],
                        help='Source model/dataset name for variable unit standardization and default input path')
    parser.add_argument('--input_dir', default=None,
                        help='Directory containing individual zarr files (default: /pscratch/sd/w/wcmca1/hackathon/etc_data/{source}/single_vars)')
    parser.add_argument('--output_dir', default=None,
                        help='Output directory (default: parent directory of input_dir)')
    parser.add_argument('--output_prefix', default='etc_2d_combined',
                        help='Prefix for output filenames (default: etc_2d_combined)')
    parser.add_argument('--suffix', default=None,
                        help='Process only files with specific suffix (e.g., track1014, all_all)')
    parser.add_argument('--chunk_size', type=int, default=1000,
                        help='Chunk size for time dimension (default: 1000)')
    parser.add_argument('--verify', action='store_true', default=True,
                        help='Verify combined files after creation (default: True)')
    parser.add_argument('--no-verify', dest='verify', action='store_false',
                        help='Skip verification step')
    parser.add_argument('--add-cof-data', action='store_true', default=True,
                        help='Add ETC COF overlap data to combined files (default: True)')
    parser.add_argument('--no-cof-data', dest='add_cof_data', action='store_false',
                        help='Skip adding COF overlap data')
    parser.add_argument('--etc-path', default=None,
                        help='Path to ETC COF parquet files (default: /pscratch/sd/w/wcmca1/hackathon/etc_tracks/)')
    
    args = parser.parse_args()
    
    # Set input directory - default based on source
    if args.input_dir is None:
        args.input_dir = f"/pscratch/sd/w/wcmca1/hackathon/etc_data/{args.source}/single_vars"
        print(f"Using default input directory for source '{args.source}'")
    
    # Verify input directory exists
    if not os.path.exists(args.input_dir):
        print(f"ERROR: Input directory does not exist: {args.input_dir}")
        return
    
    # Set output directory - default is parent directory of input_dir
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = os.path.dirname(os.path.abspath(args.input_dir))
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*70)
    print("ETC 2D Variable Combiner")
    print("="*70)
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Output prefix: {args.output_prefix}")
    if args.source:
        print(f"Source (for unit standardization): {args.source}")
    if args.add_cof_data:
        print(f"Add COF overlap data: Yes")
        if args.etc_path:
            print(f"COF data path: {args.etc_path}")
        else:
            print(f"COF data path: /pscratch/sd/w/wcmca1/hackathon/etc_tracks/ (default)")
    else:
        print(f"Add COF overlap data: No")
    if args.suffix:
        print(f"Processing suffix: {args.suffix}")
    print("="*70)
    sys.stdout.flush()
    
    # Group zarr files by suffix
    groups = group_zarr_files(args.input_dir)
    
    if not groups:
        print("ERROR: No zarr files found!")
        return
    
    print(f"\nFound {len(groups)} suffix group(s):")
    for suffix, file_list in sorted(groups.items()):
        print(f"  {suffix}: {len(file_list)} variables")
    print()
    sys.stdout.flush()
    
    # Filter by suffix if specified
    if args.suffix:
        if args.suffix in groups:
            groups = {args.suffix: groups[args.suffix]}
            print(f"Processing only suffix: {args.suffix}")
        else:
            print(f"ERROR: Suffix '{args.suffix}' not found!")
            print(f"Available suffixes: {list(groups.keys())}")
            return
    
    # Process each group
    combined_files = []
    for suffix, file_list in sorted(groups.items()):
        print(f"\n{'='*70}")
        print(f"Processing suffix group: {suffix}")
        print(f"{'='*70}")
        sys.stdout.flush()
        
        # Create output filename
        output_filename = f"{args.output_prefix}_{suffix}.zarr"
        output_path = os.path.join(output_dir, output_filename)
        
        # Combine files
        result = combine_zarr_files(file_list, output_path, suffix, 
                                   source=args.source, chunk_size=args.chunk_size,
                                   add_cof_data=args.add_cof_data, etc_path=args.etc_path)
        
        if result:
            combined_files.append(result)
            
            # Verify if requested
            if args.verify:
                verify_combined_file(result)
    
    # Print summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Successfully combined {len(combined_files)} file group(s):")
    for filepath in combined_files:
        print(f"  {filepath}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
