"""
Create ETC composites separated by hemisphere and overlap category.

This script:
1. Loads ETC COF parquet data (contains overlap flags and track IDs)
2. Loads ETC 2D Zarr data (contains environmental fields)
3. Combines the datasets by matching storm_id and time
4. Creates composites for NH and SH separately
5. Converts mask variables to binary (0/1) for frequency calculation
6. Saves composites as compressed NetCDF files

Note: Variable unit standardization is now applied in combine_etc_2d_vars.py

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import os
import sys
import argparse
import xarray as xr
import numpy as np
import pandas as pd


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Create ETC composites by hemisphere and overlap category'
    )
    parser.add_argument(
        '--source', 
        type=str, 
        required=True,
        choices=['scream', 'era5', 'nicam', 'icon', 'cesm2', 'um'],
        help='Source model/dataset name'
    )
    parser.add_argument(
        '--etc-path',
        type=str,
        default='/pscratch/sd/w/wcmca1/hackathon/etc_tracks/',
        help='Path to ETC COF parquet files'
    )
    parser.add_argument(
        '--zarr-path',
        type=str,
        default='/pscratch/sd/w/wcmca1/hackathon/etc_data/',
        help='Base path to ETC 2D Zarr files'
    )
    parser.add_argument(
        '--out-dir',
        type=str,
        default=None,
        help='Output directory (default: {zarr-path}/stats/{source}/)'
    )
    parser.add_argument(
        '--nh-lat-threshold',
        type=float,
        default=20.0,
        help='Latitude threshold for Northern Hemisphere (default: 20)'
    )
    parser.add_argument(
        '--sh-lat-threshold',
        type=float,
        default=-20.0,
        help='Latitude threshold for Southern Hemisphere (default: -20)'
    )
    
    return parser.parse_args()


def load_etc_cof_data(etc_file):
    """Load ETC COF parquet file."""
    print(f"\nLoading ETC COF data from: {etc_file}")
    if not os.path.isfile(etc_file):
        raise FileNotFoundError(f"ETC COF file not found: {etc_file}")
    
    etc_df = pd.read_parquet(etc_file)
    print(f"  Loaded {len(etc_df)} records")
    print(f"  Unique storms: {etc_df['storm_id'].nunique()}")
    
    return etc_df


def load_etc_zarr_data(zarr_file):
    """Load ETC 2D Zarr data."""
    print(f"\nLoading ETC 2D Zarr data from: {zarr_file}")
    if not os.path.isdir(zarr_file):
        raise FileNotFoundError(f"ETC Zarr directory not found: {zarr_file}")
    
    ds = xr.open_zarr(zarr_file, consolidated=True)
    print(f"  Dataset dimensions: {dict(ds.dims)}")
    print(f"  Number of variables: {len(ds.data_vars)}")
    
    return ds


def combine_datasets(ds, etc_df):
    """Combine Zarr dataset with parquet DataFrame."""
    print("\nCombining datasets...")
    
    # Ensure base_time in DataFrame matches the time coordinate in Zarr
    etc_df['time'] = pd.to_datetime(etc_df['base_time'])
    
    # Create arrays for the new variables matching Zarr's time dimension
    n_times = len(ds.time)
    
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
    
    print(f"  Created lookup dictionary with {len(etc_lookup)} entries")
    
    # Match Zarr time points with DataFrame
    for i, (time_val, storm_id_val) in enumerate(zip(ds.time.values, ds.storm_id.values)):
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
    print(f"  Matched {matched}/{n_times} time points ({100*matched/n_times:.1f}%)")
    
    print(f"  Overlap flag distribution:")
    for flag in [0, 1, 2, 3]:
        count = np.sum(overlap_flag_array == flag)
        if matched > 0:
            print(f"    Flag {flag}: {count} ({100*count/matched:.1f}% of matched)")
        else:
            print(f"    Flag {flag}: {count}")
    
    # Add the new variables to the Zarr dataset
    ds['overlap_flag'] = (['time'], overlap_flag_array, {
        'long_name': 'Co-occurrence overlap flag',
        'description': '0=isolated, 1=MCS only, 2=AR only, 3=MCS+AR',
        'units': '1'
    })
    
    # For ar_tracks and mcs_tracks, we'll store them as string representations
    ar_tracks_str = [str(tracks) if tracks is not None and len(tracks) > 0 else '[]' 
                     for tracks in ar_tracks_list]
    mcs_tracks_str = [str(tracks) if tracks is not None and len(tracks) > 0 else '[]' 
                      for tracks in mcs_tracks_list]
    
    ds['ar_tracks_str'] = (['time'], ar_tracks_str, {
        'long_name': 'AR track IDs',
        'description': 'List of AR tracks overlapping with ETC (string representation)',
        'units': '1'
    })
    
    ds['mcs_tracks_str'] = (['time'], mcs_tracks_str, {
        'long_name': 'MCS track IDs',
        'description': 'List of MCS tracks overlapping with ETC (string representation)',
        'units': '1'
    })
    
    ds['cof_lat'] = (['time'], cof_lat_array, {
        'long_name': 'ETC center latitude from COF data',
        'description': 'Latitude of ETC center from co-occurrence feature tracking',
        'units': 'degrees_north'
    })
    
    ds['cof_lon'] = (['time'], cof_lon_array, {
        'long_name': 'ETC center longitude from COF data',
        'description': 'Longitude of ETC center from co-occurrence feature tracking',
        'units': 'degrees_east'
    })
    
    print(f"  Added variables: overlap_flag, ar_tracks_str, mcs_tracks_str, cof_lat, cof_lon")
    
    return ds


def create_composites(ds, nh_lat_threshold=20.0, sh_lat_threshold=-20.0):
    """Create composites by hemisphere and overlap category."""
    print("\nCreating composites by hemisphere...")
    
    composites_nh = {}  # Northern Hemisphere
    composites_sh = {}  # Southern Hemisphere
    
    # Hemisphere masks
    mask_nh = ds.cof_lat > nh_lat_threshold
    mask_sh = ds.cof_lat < sh_lat_threshold
    
    print(f"  Northern Hemisphere (lat > {nh_lat_threshold}): {mask_nh.sum().values} points")
    print(f"  Southern Hemisphere (lat < {sh_lat_threshold}): {mask_sh.sum().values} points")
    
    # Identify mask variables that need binary conversion
    mask_vars = ['ar_mcs_etc_overlap_mask', 'mcs_ar_etc_overlap_mask', 'etc_mcs_ar_overlap_mask']
    available_mask_vars = [v for v in mask_vars if v in ds.data_vars]
    
    if available_mask_vars:
        print(f"\n  Converting mask variables to binary (0/1): {available_mask_vars}")
        
        # Create a copy of the dataset with binary masks
        ds_binary = ds.copy()
        
        # Convert masks to binary (0/1)
        ar_mask_binary = xr.where(ds['ar_mcs_etc_overlap_mask'] > 0, 1, 0) if 'ar_mcs_etc_overlap_mask' in ds else None
        mcs_mask_binary = xr.where(ds['mcs_ar_etc_overlap_mask'] > 0, 1, 0) if 'mcs_ar_etc_overlap_mask' in ds else None
        etc_mask_binary = xr.where(ds['etc_mcs_ar_overlap_mask'] > 0, 1, 0) if 'etc_mcs_ar_overlap_mask' in ds else None
        
        # Rename mask variables to frequency names
        if ar_mask_binary is not None:
            ds_binary['ar_freq'] = ar_mask_binary
            ds_binary['ar_freq'].attrs = {
                'long_name': 'AR frequency',
                'description': 'Frequency of AR occurrence (0-1)',
                'units': '1'
            }
        
        if mcs_mask_binary is not None:
            ds_binary['mcs_freq'] = mcs_mask_binary
            ds_binary['mcs_freq'].attrs = {
                'long_name': 'MCS frequency',
                'description': 'Frequency of MCS occurrence (0-1)',
                'units': '1'
            }
        
        if etc_mask_binary is not None:
            ds_binary['etc_freq'] = etc_mask_binary
            ds_binary['etc_freq'].attrs = {
                'long_name': 'ETC frequency',
                'description': 'Frequency of ETC occurrence (0-1)',
                'units': '1'
            }
        
        # Remove old mask variable names
        for old_name in available_mask_vars:
            if old_name in ds_binary:
                ds_binary = ds_binary.drop_vars(old_name)
        
        # Create exclusive precipitation variables if pr exists
        if 'pr' in ds and ar_mask_binary is not None and mcs_mask_binary is not None:
            print(f"  Creating exclusive precipitation variables (pr_ar, pr_mcs)")
            
            # Create exclusive masks (only one type present, not the other)
            # AR exclusive: AR present, but not MCS
            ar_exclusive = (ar_mask_binary == 1) & (mcs_mask_binary == 0)
            
            # MCS exclusive: MCS present, but not AR
            mcs_exclusive = (mcs_mask_binary == 1) & (ar_mask_binary == 0)
            
            # Apply masks to precipitation
            ds_binary['pr_ar'] = xr.where(ar_exclusive, ds['pr'], 0)
            ds_binary['pr_ar'].attrs = {
                'long_name': 'Precipitation under AR only',
                'description': 'Precipitation occurring under AR mask exclusively (not under MCS)',
                'units': ds['pr'].attrs.get('units', 'kg m-2 s-1')
            }
            
            ds_binary['pr_mcs'] = xr.where(mcs_exclusive, ds['pr'], 0)
            ds_binary['pr_mcs'].attrs = {
                'long_name': 'Precipitation under MCS only',
                'description': 'Precipitation occurring under MCS mask exclusively (not under AR)',
                'units': ds['pr'].attrs.get('units', 'kg m-2 s-1')
            }
    else:
        print("  No mask variables found - using original dataset")
        ds_binary = ds
    
    # Create composites for each overlap category and hemisphere
    print("\n  Creating composites for each overlap category:")
    overlap_categories = [
        ('isolated', 0, 'Isolated (no overlaps)'),
        ('mcs_only', 1, 'MCS only'),
        ('ar_only', 2, 'AR only'),
        ('3way', 3, 'MCS+AR (3-way)')
    ]
    
    for overlap_name, overlap_flag_val, description in overlap_categories:
        # NH composites
        mask_nh_overlap = (ds_binary.overlap_flag == overlap_flag_val) & mask_nh
        nh_count = mask_nh_overlap.sum().values
        composites_nh[overlap_name] = ds_binary.where(mask_nh_overlap, drop=False).mean(
            dim='time', skipna=True, keep_attrs=True
        )
        
        # SH composites
        mask_sh_overlap = (ds_binary.overlap_flag == overlap_flag_val) & mask_sh
        sh_count = mask_sh_overlap.sum().values
        composites_sh[overlap_name] = ds_binary.where(mask_sh_overlap, drop=False).mean(
            dim='time', skipna=True, keep_attrs=True
        )
        
        print(f"    {description:25s} - NH: {nh_count:6d} points, SH: {sh_count:6d} points")
    
    if available_mask_vars:
        print(f"\n  Note: Frequency variables (ar_freq, mcs_freq, etc_freq) represent occurrence frequency (0-1)")
        print(f"  Note: Exclusive precipitation variables (pr_ar, pr_mcs) are mutually exclusive")
    
    return composites_nh, composites_sh


def save_composites(composites_nh, composites_sh, out_dir):
    """Save composites as compressed NetCDF files."""
    print(f"\nSaving composites to: {out_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    
    # Compression settings
    encoding = {}
    for var in composites_nh['isolated'].data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}
    
    # Save NH composites
    print("  Saving Northern Hemisphere composites:")
    for comp_name, comp_data in composites_nh.items():
        output_path = f"{out_dir}/etc_2d_composite_nh_{comp_name}.nc"
        comp_data.to_netcdf(output_path, encoding=encoding)
        file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
        print(f"    {comp_name:10s} -> {output_path} ({file_size_mb:.1f} MB)")
    
    # Save SH composites
    print("  Saving Southern Hemisphere composites:")
    for comp_name, comp_data in composites_sh.items():
        output_path = f"{out_dir}/etc_2d_composite_sh_{comp_name}.nc"
        comp_data.to_netcdf(output_path, encoding=encoding)
        file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
        print(f"    {comp_name:10s} -> {output_path} ({file_size_mb:.1f} MB)")
    
    print(f"\n  Successfully saved {len(composites_nh) + len(composites_sh)} composite files")


def main():
    """Main function."""
    args = parse_args()
    
    print("="*80)
    print("ETC Composite Creation")
    print("="*80)
    print(f"Source: {args.source}")
    
    # Set up paths
    etc_file = f"{args.etc_path}/{args.source}_etc_cof_data.parquet"
    zarr_file = f"{args.zarr_path}/{args.source}/etc_2d_combined_all_all.zarr"
    
    if args.out_dir is None:
        out_dir = f"{args.zarr_path}/stats/{args.source}/"
    else:
        out_dir = args.out_dir
    
    try:
        # Load data
        etc_df = load_etc_cof_data(etc_file)
        ds = load_etc_zarr_data(zarr_file)
        
        # Note: Variable unit standardization is now applied in combine_etc_2d_vars.py
        # The zarr file should already have standardized units
        
        # Combine datasets
        ds = combine_datasets(ds, etc_df)
        
        # Create composites
        composites_nh, composites_sh = create_composites(
            ds, 
            nh_lat_threshold=args.nh_lat_threshold,
            sh_lat_threshold=args.sh_lat_threshold
        )
        
        # Save composites
        save_composites(composites_nh, composites_sh, out_dir)
        
        print("\n" + "="*80)
        print("SUCCESS: Composite creation completed")
        print("="*80)
        
        return 0
        
    except Exception as e:
        print(f"\nERROR: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
