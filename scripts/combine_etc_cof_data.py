#!/usr/bin/env python
"""
Combine ETC track files with COF (Co-Occurrence Feature) overlap tracking data.

This script:
1. Parses ETC track text files
2. Loads COF overlap tracking parquet files
3. Merges them by (time, lon, lat) - never by storm_id/etc_track, which is a running counter that the mask netCDF's own
   numbering need not share with this source's track file (true for ERA5: its masks are built from the full multi-decade
   tracking, so ETC_int_tag/TC_int_tag run much higher than era5's own 2019-2021 file's storm_id; see id_reference_file
   and docs/procedures/rerun_after_tracking_update.md)
4. Saves combined data to parquet format

Author: Zhe Feng
Last updated: 2026-09-28 (matched by position, not storm ID; fixes ERA5's overlap_flag coming out all-NaN)
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


def _time_key(series):
    """A Series of base_time/time as int64 nanoseconds, safe for exact-equality joins (mirrors
    extract_environments/subset_etc_env_store.py's point_keys)."""
    return pd.to_datetime(series).values.astype('datetime64[ns]').astype('int64')


def combine_etc_cof_data(etc_file, cof_file, output_file, id_reference_file=None):
    """
    Combine ETC track data with COF overlap tracking data, matched by (time, lon, lat), not by storm ID.

    Parameters:
    -----------
    etc_file : str
        Path to ETC track text file (the output's rows and columns other than overlap_flag/ar_tracks/mcs_tracks)
    cof_file : str
        Path to COF parquet file (etc_track/time/overlap_flag/ar_tracks/mcs_tracks; etc_track is the storm ID the mask
        netCDF used, which for most sources is this file's own storm_id but need not be - see id_reference_file)
    output_file : str
        Path for output combined parquet file
    id_reference_file : str, optional
        The track file whose storm_id numbering the mask (and so cof_file's etc_track) actually used, when it differs
        from etc_file - e.g. ERA5's masks are built from the full multi-decade tracking, so ETC_int_tag/TC_int_tag run
        into the tens of thousands while era5's own 2019-2021 file numbers storms from 1. Default: etc_file itself
        (true for every source except ERA5, since etc_track is the mask's own storm_id if nothing says otherwise).

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

    # Resolve each COF record's (etc_track, time) to a position, using whichever file the mask's storm IDs actually came
    # from (etc_file itself unless id_reference_file says otherwise): a small table keyed by (storm_id, time as int64 ns),
    # first occurrence wins, attached to cof_df as (_lon, _lat). Rows whose (etc_track, time) is not a point of that file
    # (should not happen for a mask built from it) get no position and so never match below.
    ref_df = etc_df if id_reference_file is None else parse_etc_track_file(id_reference_file, unstructured_mesh=unstructured_mesh)
    ref_pos = pd.DataFrame({'storm_id': ref_df['storm_id'], '_t': _time_key(ref_df['base_time']),
                             '_lon': ref_df['lon'].round(6), '_lat': ref_df['lat'].round(6)}).drop_duplicates(['storm_id', '_t'])
    cof_df = cof_df.assign(_t=_time_key(cof_df['time']))
    cof_df = cof_df.merge(ref_pos, left_on=['etc_track', '_t'], right_on=['storm_id', '_t'], how='left').drop(columns=['storm_id'])
    n_unresolved = int(cof_df['_lon'].isna().sum())
    if n_unresolved:
        print(f"  Warning: {n_unresolved} of {len(cof_df)} COF records' (etc_track, time) are not a point of "
              f"{id_reference_file or etc_file} (dropped, cannot be positioned)")
        cof_df = cof_df.dropna(subset=['_lon', '_lat'])

    # Merge by (time, lon, lat), never by storm ID (a mask's ETC_int_tag/TC_int_tag is a running counter that need not
    # match this source's own track file's storm_id - see docs/procedures/rerun_after_tracking_update.md)
    print("Merging ETC and COF data (matched by time, lon, lat, not by storm ID)...")
    etc_df = etc_df.assign(_t=_time_key(etc_df['base_time']), _lon=etc_df['lon'].round(6), _lat=etc_df['lat'].round(6))
    combined_df = etc_df.merge(
        cof_df.drop(columns=['etc_track', 'time']),
        on=['_t', '_lon', '_lat'],
        how='left'  # Keep all ETC tracks, even if no COF data
    )
    combined_clean_df = combined_df.drop(columns=['_t', '_lon', '_lat'])

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
    # Format: (etc_filename, cof_filename, output_basename, id_reference_file). id_reference_file is the track file whose
    # storm IDs the mask netCDF (and so the COF parquet's etc_track) actually used, when it is not etc_filename itself:
    # ERA5's masks are built from the full 42-year tracking (see docs/procedures/rerun_after_tracking_update.md), so its
    # ETC_int_tag/TC_int_tag numbering only matches that file, not era5's own 2019-2021 etc_filename below.
    ERA5_ID_REFERENCE = ('/global/cfs/cdirs/m1867/beharrop/kmscale_hackathon/stitch_nodes_data/'
                          'era5_tracking_full_etc_nocoldcoreonly/era5.etc_stitched_nodes.filtered_out_tcs.qs_filter_r30_d48.txt')
    file_mappings = [
        ('casesm2_10km_nocumulus_hp8.etc_stitched_nodes.filtered_out_tcs.txt',
         'casesm2_10km_nocumulus_etc_overlap_tracking.parquet',
         'casesm2_10km_nocumulus', None),

        ('era5.etc_stitched_nodes.filtered_out_tcs.txt',
         'IMERGv7_etc_overlap_tracking.parquet',
         'era5', ERA5_ID_REFERENCE),

        ('icon_d3hp003_hp8.etc_stitched_nodes.filtered_out_tcs.txt',
         'icon_d3hp003_etc_overlap_tracking.parquet',
         'icon_d3hp003', None),

        ('nicam_gl11_hp8.etc_stitched_nodes.filtered_out_tcs.txt',
         'nicam_gl11_etc_overlap_tracking.parquet',
         'nicam_gl11', None),

        ('screamv2_ne120_hp8.etc_stitched_nodes.filtered_out_tcs.txt',
         'scream_etc_overlap_tracking.parquet',
         'scream', None),

        ('um_glm_n2560_RAL3p3_hp8.etc_stitched_nodes.filtered_out_tcs.txt',
         'um_glm_n2560_RAL3p3_etc_overlap_tracking.parquet',
         'um_glm_n2560_RAL3p3', None),
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
            (etc, cof, base, ref) for etc, cof, base, ref in file_mappings
            if args.source.lower() in base.lower()
        ]
        if not file_mappings:
            print(f"Error: No files found for source '{args.source}'")
            print(f"Available sources: casesm2_10km_nocumulus, era5, icon_d3hp003, nicam_gl11, scream, um_glm_n2560_RAL3p3")
            return 1
    
    # Process each file pair
    success_count = 0
    for etc_filename, cof_filename, output_basename, id_reference_file in file_mappings:
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

        if id_reference_file and not os.path.exists(id_reference_file):
            print(f"  ERROR: ID reference file not found: {id_reference_file}")
            continue

        try:
            combined_df = combine_etc_cof_data(etc_file, cof_file, output_file, id_reference_file)
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
    # exit status 1 when an input is missing or a merge failed (they were only printed before, and the exit status was 0)
    return 0 if success_count == len(file_mappings) else 1


if __name__ == "__main__":
    sys.exit(main())
