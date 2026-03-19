#!/usr/bin/env python
"""
Combine ETC track files with COF (Co-Occurrence Feature) overlap tracking data.

This script:
1. Parses ETC track text files
2. Loads COF overlap tracking parquet files
3. Merges them based on storm_id/etc_track and base_time/time
4. Saves combined data to parquet format

Author: Zhe Feng
Last updated: November 2025
"""
import numpy as np
import pandas as pd
import os
import sys
import argparse
from pathlib import Path

# Add parent directory to path to import from src
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.env_extract_utilities import parse_etc_track_file


def combine_etc_cof_data(etc_file, cof_file, output_file):
    """
    Combine ETC track data with COF overlap tracking data.
    
    Parameters:
    -----------
    etc_file : str
        Path to ETC track text file
    cof_file : str
        Path to COF parquet file
    output_file : str
        Path for output combined parquet file
    
    Returns:
    --------
    pd.DataFrame
        Combined dataframe
    """
    # Parse ETC track file
    # ERA5 uses structured mesh (lat/lon grid), others use unstructured HEALPix
    unstructured_mesh = 'era5' not in Path(etc_file).name.lower()
    etc_df = parse_etc_track_file(etc_file, unstructured_mesh=unstructured_mesh)
    
    # Load COF overlap tracking data
    print(f"Loading COF data: {cof_file}")
    cof_df = pd.read_parquet(cof_file)
    print(f"  Loaded {len(cof_df)} COF records for {cof_df['etc_track'].nunique()} unique ETC tracks")
    
    # Merge the two dataframes
    print("Merging ETC and COF data...")
    combined_df = etc_df.merge(
        cof_df,
        left_on=['storm_id', 'base_time'],
        right_on=['etc_track', 'time'],
        how='left'  # Keep all ETC tracks, even if no COF data
    )
    
    # Clean up duplicate columns
    combined_clean_df = combined_df.drop(columns=['etc_track', 'time'])
    
    print(f"  Combined shape: {combined_clean_df.shape}")
    print(f"  ETC tracks: {combined_clean_df['storm_id'].nunique()}")
    
    # Check for missing COF data
    missing_cof = combined_clean_df['overlap_flag'].isna().sum()
    if missing_cof > 0:
        print(f"  Warning: {missing_cof} ETC points missing COF data")
    
    # Save to parquet
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    combined_clean_df.to_parquet(output_file, index=False, engine='pyarrow')
    print(f"  Saved to: {output_file}")
    
    return combined_clean_df


def main():
    parser = argparse.ArgumentParser(
        description='Combine ETC track files with COF overlap tracking data.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Process all standard files
  python combine_etc_cof_data.py
  
  # Process specific source
  python combine_etc_cof_data.py --source scream
  
  # Custom directories
  python combine_etc_cof_data.py --etc_dir /path/to/etc --cof_dir /path/to/cof --output_dir /path/to/output
        """
    )
    
    parser.add_argument('--etc_dir', default='/pscratch/sd/w/wcmca1/hackathon/etc_tracks/',
                        help='Directory containing ETC track files')
    parser.add_argument('--cof_dir', default='/pscratch/sd/w/wcmca1/hackathon/cof_masks/stats/',
                        help='Directory containing COF parquet files')
    parser.add_argument('--output_dir', default=None,
                        help='Output directory (default: same as etc_dir)')
    parser.add_argument('--source', default=None,
                        help='Process only specific source (e.g., scream, era5, nicam)')
    
    args = parser.parse_args()
    
    # Set output directory
    output_dir = args.output_dir if args.output_dir else args.etc_dir
    
    # Define file mappings
    # Format: (etc_filename, cof_filename, output_basename)
    file_mappings = [
        ('casesm2_10km_nocumulus_hp8.etc_stitched_nodes.txt', 
         'casesm2_10km_nocumulus_etc_overlap_tracking.parquet',
         'casesm2_10km_nocumulus'),
        
        ('era5.etc_stitched_nodes.txt', 
         'IMERGv7_etc_overlap_tracking.parquet',
         'era5'),
        
        ('icon_d3hp003_hp8.etc_stitched_nodes.txt', 
         'icon_d3hp003_etc_overlap_tracking.parquet',
         'icon_d3hp003'),
        
        ('nicam_gl11_hp8.etc_stitched_nodes.txt', 
         'nicam_gl11_etc_overlap_tracking.parquet',
         'nicam_gl11'),
        
        ('screamv2_ne120_hp8.etc_stitched_nodes.txt', 
         'scream_etc_overlap_tracking.parquet',
         'scream'),
        
        ('um_glm_n2560_RAL3p3_hp8.etc_stitched_nodes.txt', 
         'um_glm_n2560_RAL3p3_etc_overlap_tracking.parquet',
         'um_glm_n2560_RAL3p3'),
    ]
    
    print("="*70)
    print("ETC + COF Data Combination")
    print("="*70)
    print(f"ETC directory: {args.etc_dir}")
    print(f"COF directory: {args.cof_dir}")
    print(f"Output directory: {output_dir}")
    if args.source:
        print(f"Processing only: {args.source}")
    print("="*70)
    
    # Filter by source if specified
    if args.source:
        file_mappings = [
            (etc, cof, base) for etc, cof, base in file_mappings 
            if args.source.lower() in base.lower()
        ]
        if not file_mappings:
            print(f"Error: No files found for source '{args.source}'")
            print(f"Available sources: casesm2_10km_nocumulus, era5, icon_d3hp003, nicam_gl11, scream, um_glm_n2560_RAL3p3")
            return
    
    # Process each file pair
    success_count = 0
    for etc_filename, cof_filename, output_basename in file_mappings:
        print(f"\n{'='*70}")
        print(f"Processing: {output_basename}")
        print(f"{'='*70}")
        
        etc_file = os.path.join(args.etc_dir, etc_filename)
        cof_file = os.path.join(args.cof_dir, cof_filename)
        output_file = os.path.join(output_dir, f"{output_basename}_etc_cof_data.parquet")
        
        # Check if input files exist
        if not os.path.exists(etc_file):
            print(f"  ERROR: ETC file not found: {etc_file}")
            continue
        
        if not os.path.exists(cof_file):
            print(f"  ERROR: COF file not found: {cof_file}")
            continue
        
        try:
            combined_df = combine_etc_cof_data(etc_file, cof_file, output_file)
            success_count += 1
            print(f"  ✅ Success!")
        except Exception as e:
            print(f"  ❌ Error: {e}")
            import traceback
            traceback.print_exc()
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Successfully processed: {success_count}/{len(file_mappings)} file pairs")
    print(f"Output location: {output_dir}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
