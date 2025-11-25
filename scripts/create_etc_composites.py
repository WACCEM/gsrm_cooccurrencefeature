"""
Create ETC composites separated by hemisphere and overlap category.

This script:
1. Loads ETC 2D Zarr data (contains environmental fields and COF overlap data)
2. Creates composites for NH and SH separately
3. Converts mask variables to binary (0/1) for frequency calculation
4. Saves composites as compressed NetCDF files

Note: Variable unit standardization and COF data integration are now applied
      in combine_etc_2d_vars.py. Run that script first to prepare the zarr file.

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
        choices=['scream', 'era5', 'nicam_gl11', 'icon_d3hp003', 'casesm2_10km_nocumulus', 'um_glm_n2560_RAL3p3'],
        help='Source model/dataset name'
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
    parser.add_argument(
        '--overlap-flag',
        type=str,
        default='all',
        help='Overlap flag to process: 0=isolated, 1=mcs_only, 2=ar_only, 3=3way, or "all" (default: all)'
    )
    
    return parser.parse_args()


def load_etc_zarr_data(zarr_file):
    """Load ETC 2D Zarr data."""
    print(f"\nLoading ETC 2D Zarr data from: {zarr_file}")
    if not os.path.isdir(zarr_file):
        raise FileNotFoundError(f"ETC Zarr directory not found: {zarr_file}")
    
    ds = xr.open_zarr(zarr_file, consolidated=True)
    print(f"  Dataset dimensions: {dict(ds.sizes)}")
    print(f"  Number of variables: {len(ds.data_vars)}")
    
    return ds


def get_overlap_category(overlap_flag):
    """
    Map overlap flag to overlap category name.
    
    Parameters:
    -----------
    overlap_flag : int
        Overlap flag value (0, 1, 2, or 3)
    
    Returns:
    --------
    str : Overlap category name
    """
    overlap_map = {
        0: 'isolated',   # No overlaps - only ETC
        1: 'mcs_only',   # ETC + MCS only
        2: 'ar_only',    # ETC + AR only
        3: '3way'        # ETC + AR + MCS (3-way)
    }
    return overlap_map.get(overlap_flag, '3way')


def convert_masks_and_create_exclusive_precip(ds, overlap_category='3way'):
    """
    Convert overlap masks to binary frequency variables and create exclusive precipitation variables.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Dataset containing overlap mask variables
    overlap_category : str
        Type of overlap analysis:
        - '3way': ETC+AR+MCS (3-way overlap, overlap_flag=3)
        - 'mcs_only': ETC+MCS only (overlap_flag=1)
        - 'ar_only': ETC+AR only (overlap_flag=2)
        - 'isolated': ETC only, no overlaps (overlap_flag=0)
    
    Returns:
    --------
    ds_binary : xarray.Dataset
        Dataset with binary frequency variables and exclusive precipitation
    """
    print(f"\n  Converting masks for overlap category: {overlap_category}")
    
    ds_binary = ds.copy()
    
    if overlap_category == '3way':
        # 3-way overlap: ETC+AR+MCS
        mask_vars = ['ar_mcs_etc_overlap_mask', 'mcs_ar_etc_overlap_mask', 'etc_mcs_ar_overlap_mask']
        available_mask_vars = [v for v in mask_vars if v in ds.data_vars]
        
        if not available_mask_vars:
            print("  WARNING: No mask variables found for ETC+AR+MCS")
            return ds_binary
        
        print(f"  Converting mask variables to binary (0/1): {available_mask_vars}")
        
        # Convert masks to binary (0/1)
        ar_mask_binary = xr.where(ds['ar_mcs_etc_overlap_mask'] > 0, 1, 0) if 'ar_mcs_etc_overlap_mask' in ds else None
        mcs_mask_binary = xr.where(ds['mcs_ar_etc_overlap_mask'] > 0, 1, 0) if 'mcs_ar_etc_overlap_mask' in ds else None
        etc_mask_binary = xr.where(ds['etc_mcs_ar_overlap_mask'] > 0, 1, 0) if 'etc_mcs_ar_overlap_mask' in ds else None
        
        # Create frequency variables
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
        
        # Create exclusive precipitation variables
        if 'pr' in ds and ar_mask_binary is not None and mcs_mask_binary is not None and etc_mask_binary is not None:
            print(f"  Creating exclusive precipitation variables (pr_ar, pr_mcs, pr_etc)")
            
            # AR exclusive: AR present, but not MCS
            ar_exclusive = (ar_mask_binary == 1) & (mcs_mask_binary == 0)
            
            # MCS exclusive: MCS present, but not AR
            mcs_exclusive = (mcs_mask_binary == 1) & (ar_mask_binary == 0)
            
            # ETC exclusive: ETC present, but not in MCS nor AR
            etc_exclusive = (etc_mask_binary == 1) & (mcs_mask_binary == 0) & (ar_mask_binary == 0)
            
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
            
            ds_binary['pr_etc'] = xr.where(etc_exclusive, ds['pr'], 0)
            ds_binary['pr_etc'].attrs = {
                'long_name': 'Precipitation under ETC only',
                'description': 'Precipitation occurring under ETC mask exclusively (not under MCS nor AR)',
                'units': ds['pr'].attrs.get('units', 'kg m-2 s-1')
            }
    
    elif overlap_category == 'mcs_only':
        # 2-way overlap: ETC+MCS only
        mask_vars = ['mcs_etc_overlap_mask', 'etc_mcs_overlap_mask']
        available_mask_vars = [v for v in mask_vars if v in ds.data_vars]
        
        if not available_mask_vars:
            print("  WARNING: No mask variables found for ETC+MCS")
            return ds_binary
        
        print(f"  Converting mask variables to binary (0/1): {available_mask_vars}")
        
        # Convert masks to binary (0/1)
        mcs_mask_binary = xr.where(ds['mcs_etc_overlap_mask'] > 0, 1, 0) if 'mcs_etc_overlap_mask' in ds else None
        etc_mask_binary = xr.where(ds['etc_mcs_overlap_mask'] > 0, 1, 0) if 'etc_mcs_overlap_mask' in ds else None
        
        # Create frequency variables
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
        
        # Create precipitation variables separated by MCS mask
        if 'pr' in ds and mcs_mask_binary is not None:
            print(f"  Creating precipitation variables (pr_mcs, pr_nonmcs)")
            
            # MCS: precipitation where MCS is present
            mcs_mask = (mcs_mask_binary == 1)
            
            # Non-MCS: precipitation where MCS is not present
            nonmcs_mask = (mcs_mask_binary == 0)
            
            ds_binary['pr_mcs'] = xr.where(mcs_mask, ds['pr'], 0)
            ds_binary['pr_mcs'].attrs = {
                'long_name': 'Precipitation under MCS',
                'description': 'Precipitation occurring where MCS mask is present',
                'units': ds['pr'].attrs.get('units', 'kg m-2 s-1')
            }
            
            ds_binary['pr_nonmcs'] = xr.where(nonmcs_mask, ds['pr'], 0)
            ds_binary['pr_nonmcs'].attrs = {
                'long_name': 'Precipitation under non-MCS',
                'description': 'Precipitation occurring where MCS mask is not present',
                'units': ds['pr'].attrs.get('units', 'kg m-2 s-1')
            }
    
    elif overlap_category == 'ar_only':
        # 2-way overlap: ETC+AR only
        mask_vars = ['etc_ar_overlap_mask', 'ar_etc_overlap_mask']
        available_mask_vars = [v for v in mask_vars if v in ds.data_vars]
        
        if not available_mask_vars:
            print("  WARNING: No mask variables found for ETC+AR")
            return ds_binary
        
        print(f"  Converting mask variables to binary (0/1): {available_mask_vars}")
        
        # Convert masks to binary (0/1)
        etc_mask_binary = xr.where(ds['etc_ar_overlap_mask'] > 0, 1, 0) if 'etc_ar_overlap_mask' in ds else None
        ar_mask_binary = xr.where(ds['ar_etc_overlap_mask'] > 0, 1, 0) if 'ar_etc_overlap_mask' in ds else None
        
        # Create frequency variables
        if etc_mask_binary is not None:
            ds_binary['etc_freq'] = etc_mask_binary
            ds_binary['etc_freq'].attrs = {
                'long_name': 'ETC frequency',
                'description': 'Frequency of ETC occurrence (0-1)',
                'units': '1'
            }
        
        if ar_mask_binary is not None:
            ds_binary['ar_freq'] = ar_mask_binary
            ds_binary['ar_freq'].attrs = {
                'long_name': 'AR frequency',
                'description': 'Frequency of AR occurrence (0-1)',
                'units': '1'
            }
        
        # Remove old mask variable names
        for old_name in available_mask_vars:
            if old_name in ds_binary:
                ds_binary = ds_binary.drop_vars(old_name)
        
        # Create precipitation variables separated by AR mask
        if 'pr' in ds and ar_mask_binary is not None:
            print(f"  Creating precipitation variables (pr_ar, pr_nonar)")
            
            # AR: precipitation where AR is present
            ar_mask = (ar_mask_binary == 1)
            
            # Non-AR: precipitation where AR is not present
            nonar_mask = (ar_mask_binary == 0)
            
            ds_binary['pr_ar'] = xr.where(ar_mask, ds['pr'], 0)
            ds_binary['pr_ar'].attrs = {
                'long_name': 'Precipitation under AR',
                'description': 'Precipitation occurring where AR mask is present',
                'units': ds['pr'].attrs.get('units', 'kg m-2 s-1')
            }
            
            ds_binary['pr_nonar'] = xr.where(nonar_mask, ds['pr'], 0)
            ds_binary['pr_nonar'].attrs = {
                'long_name': 'Precipitation under non-AR',
                'description': 'Precipitation occurring where AR mask is not present',
                'units': ds['pr'].attrs.get('units', 'kg m-2 s-1')
            }
    
    elif overlap_category == 'isolated':
        # Isolated ETC: no overlaps with AR or MCS
        # For isolated cases, we only need ETC frequency
        # No exclusive precipitation variables since there are no overlaps
        print(f"  Note: 'isolated' category has ETC only, no AR/MCS overlaps")
        
        # Just copy the dataset, no mask conversion needed for isolated
        # The overlap_flag filtering will handle selecting isolated cases
        
    else:
        print(f"  WARNING: Unknown overlap_category '{overlap_category}' - using original dataset")
        return ds
    
    return ds_binary


def create_composites(ds, overlap_category='3way', overlap_flag=0, overlap_name='isolated', 
                     nh_lat_threshold=20.0, sh_lat_threshold=-20.0):
    """
    Create composites by hemisphere for a specific overlap flag.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Dataset containing environmental variables and overlap masks
    overlap_category : str
        Type of overlap analysis ('isolated', 'mcs_only', 'ar_only', or '3way')
        Determines which mask variables to create
        Should match the overlap_flag being processed
    overlap_flag : int
        Specific overlap flag to create composite for (0=isolated, 1=mcs_only, 2=ar_only, 3=3way)
    overlap_name : str
        Name for the composite output file (e.g., 'isolated', 'mcs_only', 'ar_only', '3way')
    nh_lat_threshold : float
        Latitude threshold for Northern Hemisphere (default: 20.0)
    sh_lat_threshold : float
        Latitude threshold for Southern Hemisphere (default: -20.0)
    
    Returns:
    --------
    composite_nh, composite_sh : xarray.Dataset, xarray.Dataset
        Composite datasets for NH and SH
    """
    print(f"\nCreating composite for overlap_flag={overlap_flag} ({overlap_name})...")
    
    # Hemisphere masks - compute to avoid dask array indexing issues
    mask_nh = (ds.cof_lat > nh_lat_threshold).compute()
    mask_sh = (ds.cof_lat < sh_lat_threshold).compute()
    
    print(f"  Northern Hemisphere (lat > {nh_lat_threshold}): {mask_nh.sum().values} points")
    print(f"  Southern Hemisphere (lat < {sh_lat_threshold}): {mask_sh.sum().values} points")
    
    # Convert masks and create exclusive precipitation variables
    ds_binary = convert_masks_and_create_exclusive_precip(ds, overlap_category=overlap_category)
    
    # Create composite for this specific overlap flag
    # Compute overlap masks to avoid dask array indexing issues
    mask_nh_overlap = ((ds_binary.overlap_flag == overlap_flag) & mask_nh).compute()
    nh_count = mask_nh_overlap.sum().values
    composite_nh = ds_binary.where(mask_nh_overlap, drop=True).mean(
        dim='time', skipna=True, keep_attrs=True
    )
    
    # SH composite
    mask_sh_overlap = ((ds_binary.overlap_flag == overlap_flag) & mask_sh).compute()
    sh_count = mask_sh_overlap.sum().values
    composite_sh = ds_binary.where(mask_sh_overlap, drop=True).mean(
        dim='time', skipna=True, keep_attrs=True
    )
    
    print(f"  NH: {nh_count:6d} points, SH: {sh_count:6d} points")
    
    # Print notes about frequency and exclusive precipitation variables
    freq_vars = [v for v in ds_binary.data_vars if v.endswith('_freq')]
    excl_precip_vars = [v for v in ds_binary.data_vars if v.startswith('pr_') and v != 'pr']
    
    if freq_vars:
        print(f"  Note: Frequency variables ({', '.join(freq_vars)}) represent occurrence frequency (0-1)")
    if excl_precip_vars:
        print(f"  Note: Exclusive precipitation variables ({', '.join(excl_precip_vars)}) are mutually exclusive")
    
    return composite_nh, composite_sh, overlap_name


def save_composites(composite_nh, composite_sh, overlap_name, out_dir):
    """Save composites as compressed NetCDF files."""
    print(f"\nSaving composite '{overlap_name}' to: {out_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    
    # Compression settings
    encoding = {}
    for var in composite_nh.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}
    
    # Clean up problematic attributes (like boolean types that NetCDF can't handle)
    def clean_attributes(ds):
        """Remove or convert attributes that NetCDF4 cannot handle."""
        ds_clean = ds.copy()
        removed_attrs = []
        for var_name in ds_clean.data_vars:
            attrs_to_remove = []
            for attr_name, attr_value in ds_clean[var_name].attrs.items():
                # Remove boolean attributes or convert them to int
                if isinstance(attr_value, (bool, np.bool_)):
                    attrs_to_remove.append(attr_name)
                    removed_attrs.append(f"{var_name}.{attr_name} = {attr_value} (type: {type(attr_value).__name__})")
            for attr_name in attrs_to_remove:
                del ds_clean[var_name].attrs[attr_name]
        
        if removed_attrs:
            print(f"    Removed {len(removed_attrs)} unsupported boolean attribute(s):")
            for attr_info in removed_attrs:
                print(f"      - {attr_info}")
        
        return ds_clean
    
    # Save NH composite
    output_path_nh = f"{out_dir}/etc_2d_composite_nh_{overlap_name}.nc"
    comp_nh_clean = clean_attributes(composite_nh)
    comp_nh_clean.to_netcdf(output_path_nh, encoding=encoding)
    file_size_mb_nh = os.path.getsize(output_path_nh) / (1024 * 1024)
    print(f"  NH: {output_path_nh} ({file_size_mb_nh:.1f} MB)")
    
    # Save SH composite
    output_path_sh = f"{out_dir}/etc_2d_composite_sh_{overlap_name}.nc"
    comp_sh_clean = clean_attributes(composite_sh)
    comp_sh_clean.to_netcdf(output_path_sh, encoding=encoding)
    file_size_mb_sh = os.path.getsize(output_path_sh) / (1024 * 1024)
    print(f"  SH: {output_path_sh} ({file_size_mb_sh:.1f} MB)")
    
    print(f"\n  Successfully saved 2 composite files")


def main():
    """Main function."""
    args = parse_args()
    
    print("="*80)
    print("ETC Composite Creation")
    print("="*80)
    print(f"Source: {args.source}")
    
    # Set up paths
    zarr_file = f"{args.zarr_path}/{args.source}/etc_2d_combined_all_all.zarr"
    
    if args.out_dir is None:
        out_dir = f"{args.zarr_path}/stats/{args.source}/"
    else:
        out_dir = args.out_dir
    
    try:
        # Load data
        ds = load_etc_zarr_data(zarr_file)
        
        # Note: Variable unit standardization and COF data integration
        # are now applied in combine_etc_2d_vars.py
        # The zarr file should already have standardized units and COF overlap data
        
        # Verify COF data is present
        if 'overlap_flag' not in ds:
            raise ValueError(
                "COF overlap data not found in zarr file. "
                "Please run combine_etc_2d_vars.py with --add-cof-data (default) first."
            )
        
        print(f"\nCOF overlap data found in zarr file")
        print(f"  overlap_flag: {ds.overlap_flag.dims}")
        print(f"  cof_lat: {ds.cof_lat.dims}")
        print(f"  cof_lon: {ds.cof_lon.dims}")
        
        # Define overlap flag configurations
        overlap_configs = [
            (0, 'isolated', 'Isolated (no overlaps)'),
            (1, 'mcs_only', 'MCS only'),
            (2, 'ar_only', 'AR only'),
            (3, '3way', 'MCS+AR (3-way)')
        ]
        
        # Determine which overlap flags to process
        if args.overlap_flag == 'all':
            flags_to_process = overlap_configs
            print(f"\nProcessing all overlap flags...")
        else:
            flag_val = int(args.overlap_flag)
            flags_to_process = [config for config in overlap_configs if config[0] == flag_val]
            if not flags_to_process:
                raise ValueError(f"Invalid overlap flag: {flag_val}. Must be 0, 1, 2, 3, or 'all'")
            print(f"\nProcessing overlap flag: {flag_val}")
        
        # Create composites for each requested overlap flag
        for overlap_flag, overlap_name, description in flags_to_process:
            print(f"\n{'='*60}")
            print(f"{description}")
            print(f"{'='*60}")
            
            # Determine the appropriate overlap_category for this overlap_flag
            overlap_category = get_overlap_category(overlap_flag)
            print(f"Using overlap_category: {overlap_category}")
            
            composite_nh, composite_sh, overlap_name = create_composites(
                ds, 
                overlap_category=overlap_category,
                overlap_flag=overlap_flag,
                overlap_name=overlap_name,
                nh_lat_threshold=args.nh_lat_threshold,
                sh_lat_threshold=args.sh_lat_threshold
            )
            
            # Save composites
            save_composites(composite_nh, composite_sh, overlap_name, out_dir)
        
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
