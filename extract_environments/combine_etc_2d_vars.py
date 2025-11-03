"""
Combine individual extracted ETC 2D variable zarr files into a single multi-variable zarr file.

This script:
1. Scans a directory for etc_2d_*.zarr files
2. Groups files by their suffix pattern (e.g., track1014, all_all)
3. Combines variables from the same group into a single zarr file
4. Preserves all metadata and coordinates

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


def combine_zarr_files(file_list, output_path, suffix, chunk_size=1000):
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
    chunk_size : int
        Chunk size for time dimension
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
            data_vars = [v for v in ds.data_vars if v not in 
                        ['storm_id', 'grid_id', 'storm_lat', 'storm_lon']]
            
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
    
    # Create combined dataset
    combined_ds = xr.Dataset(
        data_vars={
            **{var_name: (['time', 'y', 'x'], data.values, data.attrs) 
               for var_name, data in datasets.items()},
            'storm_id': reference_ds['storm_id'],
            'grid_id': reference_ds['grid_id'],
            'storm_lat': reference_ds['storm_lat'],
            'storm_lon': reference_ds['storm_lon']
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
    
    # Set up chunking - same as individual files
    ny = len(combined_ds['y'])
    nx = len(combined_ds['x'])
    
    encoding = {}
    for var_name in datasets.keys():
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
        
        data_vars = [v for v in ds.data_vars if v not in 
                    ['storm_id', 'grid_id', 'storm_lat', 'storm_lon']]
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
  # Combine all zarr files in a directory, grouping by suffix
  python combine_etc_2d_vars.py --input_dir /path/to/zarr/files
  
  # Specify output directory
  python combine_etc_2d_vars.py --input_dir /path/to/zarr/files --output_dir /path/to/output
  
  # Process only specific suffix pattern
  python combine_etc_2d_vars.py --input_dir /path/to/zarr/files --suffix track1014
  
  # Custom output filename prefix
  python combine_etc_2d_vars.py --input_dir /path/to/zarr/files --output_prefix etc_combined
        """
    )
    
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing individual zarr files')
    parser.add_argument('--output_dir', default=None,
                        help='Output directory (default: same as input_dir)')
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
    
    args = parser.parse_args()
    
    # Set output directory
    output_dir = args.output_dir if args.output_dir else args.input_dir
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*70)
    print("ETC 2D Variable Combiner")
    print("="*70)
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Output prefix: {args.output_prefix}")
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
        result = combine_zarr_files(file_list, output_path, suffix, args.chunk_size)
        
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
        print(f"  {os.path.basename(filepath)}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
