"""
Example script demonstrating how to read ETC overlap tracking data
from CSV or Parquet files into Pandas DataFrames.

The files contain ETC track co-occurrence information with:
- etc_track: ETC track number (integer)
- time: timestamp
- overlap_flag: 0=isolated, 1=MCS only, 2=AR only, 3=MCS+AR
- ar_tracks: AR track numbers overlapping with this ETC (list or comma-separated string)
- mcs_tracks: MCS track numbers overlapping with this ETC (list or comma-separated string)
"""

import pandas as pd
import numpy as np

# =============================================================================
# METHOD 1: Read from Parquet (RECOMMENDED - faster, preserves data types)
# =============================================================================

def read_from_parquet(parquet_file):
    """
    Read ETC tracking data from Parquet file.
    
    Advantages:
    - Much faster for large files (10-100x speedup)
    - Preserves list data types (ar_tracks and mcs_tracks are actual lists)
    - Better compression (typically 5-10x smaller file size)
    - No parsing needed
    
    Parameters
    ----------
    parquet_file : str
        Path to the parquet file
        
    Returns
    -------
    pd.DataFrame
        DataFrame with ETC overlap tracking information
    """
    print(f"Reading from Parquet: {parquet_file}")
    
    # Read parquet file
    df = pd.read_parquet(parquet_file, engine='pyarrow')
    
    print(f"  ✅ Loaded {len(df)} records")
    print(f"  Columns: {list(df.columns)}")
    print(f"  Data types: {df.dtypes.to_dict()}")
    
    # The ar_tracks and mcs_tracks are already Python lists
    print(f"\nExample record:")
    print(f"  ETC track: {df.iloc[0]['etc_track']}")
    print(f"  Time: {df.iloc[0]['time']}")
    print(f"  Overlap flag: {df.iloc[0]['overlap_flag']}")
    print(f"  AR tracks: {df.iloc[0]['ar_tracks']} (type: {type(df.iloc[0]['ar_tracks'])})")
    print(f"  MCS tracks: {df.iloc[0]['mcs_tracks']} (type: {type(df.iloc[0]['mcs_tracks'])})")
    
    return df


# =============================================================================
# METHOD 2: Read from CSV (slower, requires parsing)
# =============================================================================

def read_from_csv(csv_file):
    """
    Read ETC tracking data from CSV file.
    
    Note: ar_tracks and mcs_tracks are stored as comma-separated strings
    and need to be parsed back into lists of integers.
    
    Parameters
    ----------
    csv_file : str
        Path to the CSV file
        
    Returns
    -------
    pd.DataFrame
        DataFrame with ETC overlap tracking information
    """
    print(f"Reading from CSV: {csv_file}")
    
    # Read CSV file
    df = pd.read_csv(csv_file, index_col=0)
    
    print(f"  ✅ Loaded {len(df)} records")
    print(f"  Columns: {list(df.columns)}")
    
    # Parse comma-separated track strings back to lists of integers
    def parse_tracks(track_str):
        """Convert comma-separated string to list of integers."""
        if pd.isna(track_str) or track_str == '':
            return []
        return [int(t) for t in str(track_str).split(',')]
    
    df['ar_tracks'] = df['ar_tracks'].apply(parse_tracks)
    df['mcs_tracks'] = df['mcs_tracks'].apply(parse_tracks)
    
    print(f"\nExample record after parsing:")
    print(f"  ETC track: {df.iloc[0]['etc_track']}")
    print(f"  Time: {df.iloc[0]['time']}")
    print(f"  Overlap flag: {df.iloc[0]['overlap_flag']}")
    print(f"  AR tracks: {df.iloc[0]['ar_tracks']} (type: {type(df.iloc[0]['ar_tracks'])})")
    print(f"  MCS tracks: {df.iloc[0]['mcs_tracks']} (type: {type(df.iloc[0]['mcs_tracks'])})")
    
    return df


# =============================================================================
# EXAMPLE USAGE AND ANALYSIS
# =============================================================================

def analyze_etc_tracking(df):
    """
    Example analysis of ETC tracking data.
    
    Parameters
    ----------
    df : pd.DataFrame
        ETC tracking DataFrame
    """
    print("\n" + "="*80)
    print("ETC TRACKING ANALYSIS")
    print("="*80)
    
    # Basic statistics
    print(f"\nBasic Statistics:")
    print(f"  Total records: {len(df)}")
    print(f"  Unique ETC tracks: {df['etc_track'].nunique()}")
    print(f"  Time range: {df['time'].min()} to {df['time'].max()}")
    
    # Overlap distribution
    print(f"\nOverlap Distribution:")
    overlap_counts = df['overlap_flag'].value_counts().sort_index()
    overlap_labels = {0: 'Isolated', 1: 'MCS only', 2: 'AR only', 3: 'MCS+AR'}
    for flag, count in overlap_counts.items():
        pct = 100 * count / len(df)
        print(f"  {overlap_labels[flag]:12s} ({flag}): {count:6d} ({pct:5.1f}%)")
    
    # Count number of overlapping tracks
    df['n_ar_tracks'] = df['ar_tracks'].apply(len)
    df['n_mcs_tracks'] = df['mcs_tracks'].apply(len)
    
    print(f"\nNumber of overlapping tracks per ETC:")
    print(f"  AR tracks - mean: {df['n_ar_tracks'].mean():.2f}, max: {df['n_ar_tracks'].max()}")
    print(f"  MCS tracks - mean: {df['n_mcs_tracks'].mean():.2f}, max: {df['n_mcs_tracks'].max()}")
    
    # Find ETCs with most overlaps
    print(f"\nTop 5 ETCs with most AR overlaps:")
    top_ar = df.groupby('etc_track')['n_ar_tracks'].sum().sort_values(ascending=False).head()
    for etc_track, n_overlaps in top_ar.items():
        print(f"  ETC {etc_track}: {n_overlaps} total AR track overlaps")
    
    print(f"\nTop 5 ETCs with most MCS overlaps:")
    top_mcs = df.groupby('etc_track')['n_mcs_tracks'].sum().sort_values(ascending=False).head()
    for etc_track, n_overlaps in top_mcs.items():
        print(f"  ETC {etc_track}: {n_overlaps} total MCS track overlaps")
    
    # Example: Filter for specific overlap type
    print(f"\nExample Queries:")
    
    # ETCs with both MCS and AR overlap
    both_overlap = df[df['overlap_flag'] == 3]
    print(f"  Records with both MCS+AR overlap: {len(both_overlap)}")
    
    # ETCs overlapping with specific AR track (e.g., AR track 100)
    ar_track_id = 100
    if any(df['ar_tracks'].apply(lambda x: ar_track_id in x)):
        with_ar_100 = df[df['ar_tracks'].apply(lambda x: ar_track_id in x)]
        print(f"  Records where ETC overlaps with AR track {ar_track_id}: {len(with_ar_100)}")
    
    # ETCs with multiple MCS overlaps at same time
    multi_mcs = df[df['n_mcs_tracks'] > 1]
    print(f"  Records with multiple MCS overlaps: {len(multi_mcs)}")


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """Main function demonstrating file reading."""
    
    # Example file paths - UPDATE THESE TO YOUR ACTUAL FILES
    source_name = "scream"  # Change this to your source name
    output_dir = "/pscratch/sd/w/wcmca1/hackathon/cof_masks/test/"
    
    parquet_file = f"{output_dir}/{source_name}_etc_coftracks_hp8_v1.parquet"
    csv_file = f"{output_dir}/{source_name}_etc_coftracks_hp8_v1.csv"
    
    print("="*80)
    print("ETC TRACKING DATA READER EXAMPLE")
    print("="*80)
    
    # Try reading from Parquet first (recommended)
    try:
        df = read_from_parquet(parquet_file)
        analyze_etc_tracking(df)
        
    except FileNotFoundError:
        print(f"\n⚠️  Parquet file not found: {parquet_file}")
        print("Trying CSV instead...")
        
        try:
            df = read_from_csv(csv_file)
            analyze_etc_tracking(df)
            
        except FileNotFoundError:
            print(f"\n❌ CSV file not found either: {csv_file}")
            print("\nPlease update the file paths in the script to match your files.")
            return
    
    # Additional example: Convert time to datetime for time series analysis
    print("\n" + "="*80)
    print("TIME SERIES EXAMPLE")
    print("="*80)
    
    # Convert time to datetime
    df['datetime'] = pd.to_datetime(df['time'])
    df = df.set_index('datetime')
    
    # Count overlaps by day
    daily_overlaps = df.groupby([df.index.date, 'overlap_flag']).size().unstack(fill_value=0)
    print("\nDaily overlap counts (first 10 days):")
    print(daily_overlaps.head(10))
    
    # Get all records for a specific ETC track
    example_etc = df['etc_track'].iloc[0]
    etc_timeseries = df[df['etc_track'] == example_etc].sort_index()
    print(f"\nTime series for ETC track {example_etc} ({len(etc_timeseries)} timesteps):")
    print(etc_timeseries[['etc_track', 'overlap_flag', 'ar_tracks', 'mcs_tracks']].head(10))


if __name__ == "__main__":
    main()
