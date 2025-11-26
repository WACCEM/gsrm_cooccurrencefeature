#!/usr/bin/env python
"""
Calculate spatial statistics for ETC 2D data.

This script computes statistics over spatial dimensions (x, y) to condense 3D data
(time, y, x) into 1D time series. Statistics include:
- Mean/min/max for regular variables within circular radius
- Fractional area coverage for feature masks
- Domain-mean precipitation separated by features (AR, MCS, ETC)

The script processes all time points uniformly, preserving overlap_flag as a coordinate
for filtering during analysis.

Author: Zhe Feng (zhe.feng@pnnl.gov)
Date: 2025-11-25
"""

import os
import sys
import argparse
import time
import xarray as xr
import numpy as np


def create_circular_mask(x_dim, y_dim, radius_deg):
    """
    Create a boolean circular mask centered at (0, 0).
    
    Parameters:
    -----------
    x_dim : xarray.DataArray
        X coordinate array (centered at 0)
    y_dim : xarray.DataArray
        Y coordinate array (centered at 0)
    radius_deg : float
        Radius in degrees
    
    Returns:
    --------
    mask : xarray.DataArray
        Boolean mask (True inside circle, False outside)
    """
    # Create distance field from center
    distance = np.sqrt(x_dim**2 + y_dim**2)
    mask = distance <= radius_deg
    return mask


def compute_spatial_stats_basic(ds, variables, radius_deg=10):
    """
    Compute mean, min, max for regular variables within circular radius.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Input dataset with dimensions (time, y, x)
    variables : list of str
        Variable names to compute statistics for
    radius_deg : float
        Radius in degrees (default: 10)
    
    Returns:
    --------
    stats_ds : xarray.Dataset
        Dataset with computed statistics (time dimension only)
    """
    # Create circular mask
    mask = create_circular_mask(ds.x, ds.y, radius_deg)
    
    # Initialize result dictionary
    stats_dict = {}
    
    for var_name in variables:
        if var_name not in ds:
            print(f"  WARNING: Variable '{var_name}' not found, skipping")
            continue
        
        # Apply mask and compute statistics
        var_masked = ds[var_name].where(mask)
        
        # Compute statistics along spatial dimensions (vectorized)
        stats_dict[f'{var_name}_mean'] = var_masked.mean(dim=['x', 'y'], skipna=True)
        stats_dict[f'{var_name}_min'] = var_masked.min(dim=['x', 'y'], skipna=True)
        stats_dict[f'{var_name}_max'] = var_masked.max(dim=['x', 'y'], skipna=True)
        
        # Add attributes
        orig_units = ds[var_name].attrs.get('units', '1')
        stats_dict[f'{var_name}_mean'].attrs = {
            'long_name': f'Mean {var_name} within {radius_deg}° radius',
            'units': orig_units
        }
        stats_dict[f'{var_name}_min'].attrs = {
            'long_name': f'Min {var_name} within {radius_deg}° radius',
            'units': orig_units
        }
        stats_dict[f'{var_name}_max'].attrs = {
            'long_name': f'Max {var_name} within {radius_deg}° radius',
            'units': orig_units
        }
    
    # Create output dataset
    stats_ds = xr.Dataset(stats_dict)
    
    return stats_ds


def compute_mask_fractional_area(ds, mask_variables, radii=[10, 15]):
    """
    Compute fractional area coverage for feature masks within circular radii.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Input dataset with dimensions (time, y, x)
    mask_variables : list of str
        Variable names ending with '_overlap_mask'
    radii : list of float
        List of radii in degrees (default: [10, 15])
    
    Returns:
    --------
    stats_ds : xarray.Dataset
        Dataset with fractional area coverage (time dimension only)
    """
    stats_dict = {}
    
    for radius_deg in radii:
        # Create circular mask
        circle_mask = create_circular_mask(ds.x, ds.y, radius_deg)
        total_pixels = circle_mask.sum().values  # Total pixels in circle
        
        for var_name in mask_variables:
            if var_name not in ds:
                print(f"  WARNING: Variable '{var_name}' not found, skipping")
                continue
            
            # Count pixels where mask > 0 within the circle
            # Vectorized: (mask > 0) & circle_mask, then sum over x,y
            feature_present = (ds[var_name] > 0) & circle_mask
            feature_count = feature_present.sum(dim=['x', 'y'])
            
            # Compute fractional area
            fractional_area = feature_count / total_pixels
            
            # Store result
            output_name = f'{var_name}_frac_r{int(radius_deg)}'
            stats_dict[output_name] = fractional_area
            stats_dict[output_name].attrs = {
                'long_name': f'Fractional area coverage of {var_name} within {radius_deg}° radius',
                'units': '1',
                'description': f'Fraction of pixels with mask > 0 within {radius_deg}° radius'
            }
    
    # Create output dataset
    stats_ds = xr.Dataset(stats_dict)
    
    return stats_ds


def compute_feature_precipitation_stats(ds, radius_deg=10):
    """
    Compute domain-mean precipitation separated by features uniformly.
    
    This computes feature-specific precipitation statistics regardless of overlap_flag.
    Variables that don't exist will naturally result in NaN or 0, which is informative.
    
    Note: Total pr_mean should be computed via compute_spatial_stats_basic() to avoid duplication.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Input dataset with 'pr' and mask variables
    radius_deg : float
        Radius in degrees (default: 10)
    
    Returns:
    --------
    stats_ds : xarray.Dataset
        Dataset with feature-specific precipitation statistics
    """
    # Create circular mask
    circle_mask = create_circular_mask(ds.x, ds.y, radius_deg)
    total_pixels = circle_mask.sum().values  # Total pixels in circle
    
    stats_dict = {}
    
    if 'pr' not in ds:
        print("  WARNING: 'pr' variable not found")
        return xr.Dataset()
    
    # Get precipitation variable
    pr = ds['pr'].where(circle_mask)  # Only consider pixels within radius
    
    # AR-related precipitation (if AR mask exists)
    if 'ar_etc_overlap_mask' in ds:
        ar_mask = (ds['ar_etc_overlap_mask'] > 0) & circle_mask
        
        # Precipitation where AR is present
        pr_ar = (pr * ar_mask).sum(dim=['x', 'y']) / total_pixels
        stats_dict['pr_ar_mean'] = pr_ar
        stats_dict['pr_ar_mean'].attrs = {
            'long_name': 'Domain-mean precipitation from AR',
            'units': pr.attrs.get('units', 'mm h-1'),
            'description': f'Sum of pr under AR mask / total pixels in {radius_deg}° radius'
        }
        
        # Precipitation where AR is not present
        pr_nonar = (pr * (~ar_mask & circle_mask)).sum(dim=['x', 'y']) / total_pixels
        stats_dict['pr_nonar_mean'] = pr_nonar
        stats_dict['pr_nonar_mean'].attrs = {
            'long_name': 'Domain-mean precipitation from non-AR',
            'units': pr.attrs.get('units', 'mm h-1'),
            'description': f'Sum of pr outside AR mask / total pixels in {radius_deg}° radius'
        }
    
    # Alternative AR mask (for 3-way overlaps)
    if 'ar_mcs_etc_overlap_mask' in ds:
        ar_mask_3way = (ds['ar_mcs_etc_overlap_mask'] > 0) & circle_mask
        
        # AR exclusive (not under MCS) - for 3-way comparisons
        if 'mcs_ar_etc_overlap_mask' in ds:
            mcs_mask_3way = (ds['mcs_ar_etc_overlap_mask'] > 0) & circle_mask
            ar_exclusive = ar_mask_3way & ~mcs_mask_3way
            
            pr_ar_excl = (pr * ar_exclusive).sum(dim=['x', 'y']) / total_pixels
            stats_dict['pr_ar_excl_mean'] = pr_ar_excl
            stats_dict['pr_ar_excl_mean'].attrs = {
                'long_name': 'Domain-mean precipitation from AR (exclusive of MCS)',
                'units': pr.attrs.get('units', 'mm h-1'),
                'description': f'Sum of pr under AR mask (not MCS) / total pixels in {radius_deg}° radius'
            }
    
    # MCS-related precipitation (if MCS mask exists)
    if 'mcs_etc_overlap_mask' in ds:
        mcs_mask = (ds['mcs_etc_overlap_mask'] > 0) & circle_mask
        
        # Precipitation where MCS is present
        pr_mcs = (pr * mcs_mask).sum(dim=['x', 'y']) / total_pixels
        stats_dict['pr_mcs_mean'] = pr_mcs
        stats_dict['pr_mcs_mean'].attrs = {
            'long_name': 'Domain-mean precipitation from MCS',
            'units': pr.attrs.get('units', 'mm h-1'),
            'description': f'Sum of pr under MCS mask / total pixels in {radius_deg}° radius'
        }
        
        # Precipitation where MCS is not present
        pr_nonmcs = (pr * (~mcs_mask & circle_mask)).sum(dim=['x', 'y']) / total_pixels
        stats_dict['pr_nonmcs_mean'] = pr_nonmcs
        stats_dict['pr_nonmcs_mean'].attrs = {
            'long_name': 'Domain-mean precipitation from non-MCS',
            'units': pr.attrs.get('units', 'mm h-1'),
            'description': f'Sum of pr outside MCS mask / total pixels in {radius_deg}° radius'
        }
    
    # Alternative MCS mask (for 3-way overlaps)
    if 'mcs_ar_etc_overlap_mask' in ds:
        mcs_mask_3way = (ds['mcs_ar_etc_overlap_mask'] > 0) & circle_mask
        
        # MCS exclusive (not under AR) - for 3-way comparisons
        if 'ar_mcs_etc_overlap_mask' in ds:
            ar_mask_3way = (ds['ar_mcs_etc_overlap_mask'] > 0) & circle_mask
            mcs_exclusive = mcs_mask_3way & ~ar_mask_3way
            
            pr_mcs_excl = (pr * mcs_exclusive).sum(dim=['x', 'y']) / total_pixels
            stats_dict['pr_mcs_excl_mean'] = pr_mcs_excl
            stats_dict['pr_mcs_excl_mean'].attrs = {
                'long_name': 'Domain-mean precipitation from MCS (exclusive of AR)',
                'units': pr.attrs.get('units', 'mm h-1'),
                'description': f'Sum of pr under MCS mask (not AR) / total pixels in {radius_deg}° radius'
            }
    
    # ETC-related precipitation (for 3-way overlaps)
    if 'etc_mcs_ar_overlap_mask' in ds:
        etc_mask = (ds['etc_mcs_ar_overlap_mask'] > 0) & circle_mask
        
        # ETC exclusive (not under MCS nor AR) - for 3-way comparisons
        if 'ar_mcs_etc_overlap_mask' in ds and 'mcs_ar_etc_overlap_mask' in ds:
            ar_mask_3way = (ds['ar_mcs_etc_overlap_mask'] > 0) & circle_mask
            mcs_mask_3way = (ds['mcs_ar_etc_overlap_mask'] > 0) & circle_mask
            etc_exclusive = etc_mask & ~mcs_mask_3way & ~ar_mask_3way
            
            pr_etc = (pr * etc_exclusive).sum(dim=['x', 'y']) / total_pixels
            stats_dict['pr_etc_mean'] = pr_etc
            stats_dict['pr_etc_mean'].attrs = {
                'long_name': 'Domain-mean precipitation from ETC (exclusive of MCS and AR)',
                'units': pr.attrs.get('units', 'mm h-1'),
                'description': f'Sum of pr under ETC mask (not MCS/AR) / total pixels in {radius_deg}° radius'
            }
    
    # Create output dataset
    stats_ds = xr.Dataset(stats_dict)
    
    return stats_ds


def compute_all_spatial_stats(ds, 
                               basic_vars=None, 
                               mask_vars=None,
                               basic_radius=10,
                               mask_radii=[10, 15],
                               pr_radius=10):
    """
    Main function to compute all spatial statistics uniformly for all time points.
    
    This function computes statistics for all time points regardless of overlap_flag.
    The overlap_flag is preserved as a coordinate variable for filtering during analysis.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Input dataset with dimensions (time, y, x)
    basic_vars : list of str, optional
        Variables for basic statistics (mean/min/max)
    mask_vars : list of str, optional
        Mask variables for fractional area
    basic_radius : float
        Radius for basic statistics (default: 10)
    mask_radii : list of float
        Radii for mask statistics (default: [10, 15])
    pr_radius : float
        Radius for precipitation statistics (default: 10)
    
    Returns:
    --------
    combined_ds : xarray.Dataset
        Combined dataset with all statistics and overlap_flag as coordinate
    """
    print(f"\nComputing spatial statistics for all time points")
    print(f"  Input shape: {dict(ds.sizes)}")
    
    stats_datasets = []
    
    # 1. Basic statistics for regular variables
    if basic_vars:
        print(f"  Computing basic stats (mean/min/max) for {len(basic_vars)} variables...")
        stats_basic = compute_spatial_stats_basic(ds, basic_vars, radius_deg=basic_radius)
        stats_datasets.append(stats_basic)
    
    # 2. Fractional area for mask variables
    if mask_vars:
        print(f"  Computing fractional area for {len(mask_vars)} mask variables...")
        stats_masks = compute_mask_fractional_area(ds, mask_vars, radii=mask_radii)
        stats_datasets.append(stats_masks)
    
    # 3. Feature precipitation statistics (computed uniformly for all time points)
    print(f"  Computing feature-specific precipitation statistics...")
    stats_pr = compute_feature_precipitation_stats(ds, radius_deg=pr_radius)
    stats_datasets.append(stats_pr)
    
    # Combine all statistics
    combined_ds = xr.merge(stats_datasets)
    
    # Add overlap_flag as a coordinate variable for easy filtering
    if 'overlap_flag' in ds:
        combined_ds = combined_ds.assign_coords({'overlap_flag': ds.overlap_flag})
        print(f"  Added overlap_flag as coordinate for filtering")
    
    # Add other useful coordinates
    if 'cof_lat' in ds:
        combined_ds = combined_ds.assign_coords({'cof_lat': ds.cof_lat})
    if 'cof_lon' in ds:
        combined_ds = combined_ds.assign_coords({'cof_lon': ds.cof_lon})
    if 'storm_id' in ds:
        combined_ds = combined_ds.assign_coords({'storm_id': ds.storm_id})
    
    # Add metadata
    combined_ds.attrs['basic_radius_deg'] = basic_radius
    combined_ds.attrs['mask_radii_deg'] = mask_radii
    combined_ds.attrs['pr_radius_deg'] = pr_radius
    combined_ds.attrs['description'] = 'Spatial statistics computed within circular radii centered at ETC. overlap_flag coordinate: 0=isolated, 1=mcs_only, 2=ar_only, 3=3way'
    
    print(f"  Output variables: {len(combined_ds.data_vars)}")
    print(f"  Output shape: {dict(combined_ds.sizes)}")
    print(f"  Coordinates: {list(combined_ds.coords)}")
    
    return combined_ds


def main():
    """Main function to process ETC spatial statistics."""
    
    parser = argparse.ArgumentParser(
        description='Calculate spatial statistics for ETC 2D data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--source', type=str, required=True,
                        help='Data source name (e.g., era5, scream, etc.)')
    parser.add_argument('--zarr-path', type=str, 
                        default='/pscratch/sd/w/wcmca1/hackathon/etc_data',
                        help='Base path to zarr files')
    parser.add_argument('--output-dir', type=str,
                        default='/pscratch/sd/w/wcmca1/hackathon/etc_data/stats',
                        help='Output directory for statistics')
    parser.add_argument('--basic-radius', type=float, default=10.0,
                        help='Radius (degrees) for basic statistics')
    parser.add_argument('--mask-radii', type=float, nargs='+', default=[10.0, 15.0],
                        help='Radii (degrees) for mask fractional area')
    parser.add_argument('--pr-radius', type=float, default=10.0,
                        help='Radius (degrees) for precipitation statistics')
    
    args = parser.parse_args()
    
    print("=" * 80)
    print(f"ETC Spatial Statistics Calculator")
    print("=" * 80)
    print(f"Source: {args.source}")
    print(f"Zarr path: {args.zarr_path}")
    print(f"Output directory: {args.output_dir}")
    print(f"Basic stats radius: {args.basic_radius}°")
    print(f"Mask radii: {args.mask_radii}°")
    print(f"Precipitation radius: {args.pr_radius}°")
    print("=" * 80)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Construct zarr file path
    zarr_file = os.path.join(args.zarr_path, args.source, 'etc_2d_combined_all_all.zarr')
    
    if not os.path.isdir(zarr_file):
        print(f"ERROR: Zarr file not found: {zarr_file}")
        sys.exit(1)
    
    print(f"\nLoading zarr file: {zarr_file}")
    start_time = time.time()
    
    # Load dataset
    ds = xr.open_zarr(zarr_file, consolidated=True)
    print(f"  Loaded in {time.time() - start_time:.2f} seconds")
    print(f"  Dataset shape: {dict(ds.sizes)}")
    print(f"  Number of variables: {len(ds.data_vars)}")
    
    # Identify variable types
    metadata_vars = ['storm_id', 'grid_id', 'lon_id', 'lat_id', 'storm_lat', 'storm_lon', 
                     'overlap_flag', 'ar_tracks_str', 'mcs_tracks_str', 'cof_lat', 'cof_lon']
    mask_vars = [v for v in ds.data_vars if v.endswith('_overlap_mask')]
    basic_vars = [v for v in ds.data_vars if v not in metadata_vars and v not in mask_vars]
    
    print(f"\n  Basic variables ({len(basic_vars)}): {', '.join(basic_vars[:5])}{'...' if len(basic_vars) > 5 else ''}")
    print(f"  Mask variables ({len(mask_vars)}): {', '.join(mask_vars)}")
    
    # Compute statistics
    print("\n" + "=" * 80)
    start_time = time.time()
    
    stats = compute_all_spatial_stats(
        ds,
        basic_vars=basic_vars,
        mask_vars=mask_vars,
        basic_radius=args.basic_radius,
        mask_radii=args.mask_radii,
        pr_radius=args.pr_radius
    )
    
    compute_time = time.time() - start_time
    print(f"\n  Statistics computed in {compute_time:.2f} seconds")
    
    # Save to NetCDF
    output_file = os.path.join(args.output_dir, f'etc_spatial_stats_{args.source}.nc')
    print(f"\nSaving to: {output_file}")
    
    start_time = time.time()
    stats.to_netcdf(output_file)
    save_time = time.time() - start_time
    
    file_size_mb = os.path.getsize(output_file) / 1e6
    print(f"  Saved in {save_time:.2f} seconds")
    print(f"  File size: {file_size_mb:.2f} MB")
    
    # Calculate compression ratio
    original_size_mb = ds.nbytes / 1e6
    compression_ratio = original_size_mb / file_size_mb
    print(f"  Original data size: {original_size_mb:.2f} MB")
    print(f"  Compression ratio: {compression_ratio:.2f}x")
    
    print("\n" + "=" * 80)
    print("Processing complete!")
    print("=" * 80)
    print("\nTo filter during analysis, use:")
    print(f"  import xarray as xr")
    print(f"  stats = xr.open_dataset('{output_file}')")
    print(f"  stats_ar_only = stats.where(stats.overlap_flag == 2, drop=True)")
    print(f"  stats_3way = stats.where(stats.overlap_flag == 3, drop=True)")
    print("=" * 80)


if __name__ == '__main__':
    main()
