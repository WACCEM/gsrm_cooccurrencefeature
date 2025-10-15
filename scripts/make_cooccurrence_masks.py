#!/usr/bin/env python3
"""
Make Co-occurrence Feature Masks Processing Script

This script processes overlaps between MCS, AR, and ETC features across multiple time steps,
defines co-occurrence features, and saves results to zarr format.

Key Features:
- Define co-occurrence features based on atmospheric feature scale hierarchy
- Consider 2-way & 3-way overlaps
- Vectorized overlap detection for efficiency
- Processing of full time series datasets using streaming to minimize memory usage
- Output in zarr format with proper coordinates and attributes

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import os
import sys
import numpy as np
import xarray as xr
from pathlib import Path
import warnings
import argparse
import logging
# import easygems.healpix as egh  # Commented out for testing
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from src.zarr_tools import stream_process_to_zarr, initialize_zarr_store, setup_dask_client

warnings.filterwarnings('ignore')

def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def find_overlapping_tracks_and_pairs(mask1, mask2, binary_sum_mask, 
                                     thresh1=0.1, thresh2=0.1, overlap_threshold=1,
                                     feature1_name="Feature1", feature2_name="Feature2", verbose=True):
    """
    Find overlapping tracks for both feature types AND their valid pairs in a single pass.
    
    Parameters:
    -----------
    mask1, mask2 : xarray.DataArray
        Track masks for the two feature types to analyze
    binary_sum_mask : xarray.DataArray
        Binary summation mask indicating overlap regions
    thresh1, thresh2 : float
        Minimum overlap fractions for each feature type
    overlap_threshold : int, default=1
        Minimum value in binary_sum_mask for overlap detection
    feature1_name, feature2_name : str
        Names of the feature types for reporting
    verbose : bool, default=True
        Whether to print progress information
        
    Returns:
    --------
    dict : Combined results containing:
        - 'feature1_results': overlap results for feature1 (like find_overlapping_tracks)
        - 'feature2_results': overlap results for feature2 (like find_overlapping_tracks) 
        - 'valid_pairs': list of (track1_id, track2_id) tuples
        - 'feature1_tracks': numpy array of feature1 track IDs in valid pairs
        - 'feature2_tracks': numpy array of feature2 track IDs in valid pairs
        - 'pair_details': detailed information about each pair
    """
    if verbose:
        print(f"Finding {feature1_name}-{feature2_name} overlaps and pairs with thresholds: {thresh1:.1%}, {thresh2:.1%}")
    
    # Step 1: Find overlap condition (shared between both features)
    overlap_condition = (binary_sum_mask > overlap_threshold) & ~np.isnan(binary_sum_mask) & ~np.isnan(mask1) & ~np.isnan(mask2)
    
    # Step 2: Get overlapping tracks for both features in one pass
    # Feature 1 overlapping tracks
    tracks1_in_overlap, overlap1_counts = np.unique(mask1.where(overlap_condition), return_counts=True)
    valid_mask1 = ~np.isnan(tracks1_in_overlap)
    tracks1_in_overlap = tracks1_in_overlap[valid_mask1]
    overlap1_counts = overlap1_counts[valid_mask1]
    
    # Feature 2 overlapping tracks
    tracks2_in_overlap, overlap2_counts = np.unique(mask2.where(overlap_condition), return_counts=True)
    valid_mask2 = ~np.isnan(tracks2_in_overlap)
    tracks2_in_overlap = tracks2_in_overlap[valid_mask2]
    overlap2_counts = overlap2_counts[valid_mask2]
    
    # Step 3: Get total pixel counts for all tracks (shared computation)
    unique_tracks1, total1_counts = np.unique(mask1, return_counts=True)
    unique_tracks2, total2_counts = np.unique(mask2, return_counts=True)
    
    track1_counts_dict = dict(zip(unique_tracks1, total1_counts))
    track2_counts_dict = dict(zip(unique_tracks2, total2_counts))
    
    # Step 4: Calculate overlap fractions for both features
    total1_pixel_counts = np.array([track1_counts_dict.get(track, 0) for track in tracks1_in_overlap])
    total2_pixel_counts = np.array([track2_counts_dict.get(track, 0) for track in tracks2_in_overlap])
    
    overlap1_fractions = overlap1_counts / total1_pixel_counts
    overlap2_fractions = overlap2_counts / total2_pixel_counts
    
    # Step 5: Find tracks exceeding overlap thresholds
    threshold1_mask = overlap1_fractions >= thresh1
    threshold2_mask = overlap2_fractions >= thresh2
    
    tracks1_exceeding = tracks1_in_overlap[threshold1_mask]
    tracks2_exceeding = tracks2_in_overlap[threshold2_mask]
    
    # Step 6: EFFICIENT PAIR FINDING using vectorized operations
    # Create masks for tracks that exceed thresholds
    mask1_qualifying = np.isin(mask1, tracks1_exceeding)
    mask2_qualifying = np.isin(mask2, tracks2_exceeding)
    
    # Find pixels where both qualifying tracks overlap
    mutual_overlap_pixels = mask1_qualifying & mask2_qualifying & overlap_condition
    
    if mutual_overlap_pixels.sum().compute().item() > 0:
        # Get track pairs that actually overlap spatially
        track1_values = mask1.where(mutual_overlap_pixels).values
        track2_values = mask2.where(mutual_overlap_pixels).values
        
        # Remove NaN values and create coordinate pairs
        valid_overlap = ~np.isnan(track1_values) & ~np.isnan(track2_values)
        track1_valid = track1_values[valid_overlap]
        track2_valid = track2_values[valid_overlap]
        
        # Find unique pairs using numpy operations
        pair_coords = np.column_stack([track1_valid, track2_valid])
        unique_pairs = np.unique(pair_coords, axis=0)
        
        valid_pairs = [(int(pair[0]), int(pair[1])) for pair in unique_pairs]
    else:
        valid_pairs = []
    
    # Step 7: Create detailed results (similar to original functions)
    overlap1_fractions_dict = {tracks1_in_overlap[i]: overlap1_fractions[i] for i in range(len(tracks1_in_overlap))}
    overlap2_fractions_dict = {tracks2_in_overlap[i]: overlap2_fractions[i] for i in range(len(tracks2_in_overlap))}
    
    # Create pair details
    pair_details = {}
    for track1, track2 in valid_pairs:
        # Count overlap pixels for this specific pair
        track1_pixels = (mask1 == track1)
        track2_pixels = (mask2 == track2)
        pair_overlap = track1_pixels & track2_pixels & overlap_condition
        overlap_pixel_count = pair_overlap.sum().compute().item()
        
        pair_details[(track1, track2)] = {
            'overlap_pixels': overlap_pixel_count,
            'track1_fraction': overlap1_fractions_dict.get(track1, 0),
            'track2_fraction': overlap2_fractions_dict.get(track2, 0)
        }
    
    # Extract unique track lists from valid pairs
    if valid_pairs:
        paired_feature1_tracks = np.unique([pair[0] for pair in valid_pairs])
        paired_feature2_tracks = np.unique([pair[1] for pair in valid_pairs])
    else:
        paired_feature1_tracks = np.array([])
        paired_feature2_tracks = np.array([])
    
    if verbose:
        print(f"Feature1 ({feature1_name}): {len(tracks1_exceeding)} tracks exceed {thresh1:.1%} threshold")
        print(f"Feature2 ({feature2_name}): {len(tracks2_exceeding)} tracks exceed {thresh2:.1%} threshold")
        print(f"Found {len(valid_pairs)} valid pairs with spatial overlap")
        if len(valid_pairs) > 0 and len(valid_pairs) <= 10:
            print(f"  Pairs: {valid_pairs}")
    
    return {
        'feature1_results': {
            'tracks_exceeding_threshold': tracks1_exceeding,
            'overlap_fractions': overlap1_fractions_dict,
            'total_tracks_analyzed': len(unique_tracks1),
            'threshold_info': {'min_overlap_fraction': thresh1, 'overlap_threshold': overlap_threshold}
        },
        'feature2_results': {
            'tracks_exceeding_threshold': tracks2_exceeding,
            'overlap_fractions': overlap2_fractions_dict,
            'total_tracks_analyzed': len(unique_tracks2),
            'threshold_info': {'min_overlap_fraction': thresh2, 'overlap_threshold': overlap_threshold}
        },
        'valid_pairs': valid_pairs,
        'feature1_tracks': paired_feature1_tracks,
        'feature2_tracks': paired_feature2_tracks,
        'pair_details': pair_details
    }


def find_true_3way_overlaps(mcs_mask, ar_mask, etc_mask, binary_sum_mask,
                          mcs_thresh=0.2, ar_thresh=0.1, etc_thresh=0.0, verbose=True):
    """
    Find tracks that have true 3-way overlaps with mutual spatial verification.
    
    This function ensures that all three features actually overlap spatially,
    not just that they individually exceed thresholds in 3-way regions.
    
    Parameters:
    -----------
    mcs_mask, ar_mask, etc_mask : xarray.DataArray
        Track masks for the three feature types
    binary_sum_mask : xarray.DataArray
        Binary summation mask where values == 3 indicate true 3-way overlaps
    mcs_thresh, ar_thresh, etc_thresh : float
        Minimum overlap fractions for each feature type
    verbose : bool, default=True
        Whether to print progress information
        
    Returns:
    --------
    dict : Results containing validated 3-way overlaps with mutual verification
    """
    if verbose:
        print("Finding true 3-way overlaps with mutual spatial verification...")
    
    # Create 3-way overlap mask (only where all 3 features overlap)
    threeway_mask = (binary_sum_mask == 3).astype(int)
    threeway_pixels = threeway_mask.sum().compute().item()
    
    if verbose:
        print(f"Total 3-way overlap pixels: {threeway_pixels}")
    
    if threeway_pixels == 0:
        return {
            'mcs_tracks': np.array([]),
            'ar_tracks': np.array([]),
            'etc_tracks': np.array([]),
            'validated_triplets': [],
            'triplet_details': {}
        }
    
    # Find all tracks that have any presence in 3-way regions
    mcs_in_3way = np.unique(mcs_mask.values[threeway_mask > 0])
    mcs_in_3way = mcs_in_3way[mcs_in_3way > 0]
    
    ar_in_3way = np.unique(ar_mask.values[threeway_mask > 0])
    ar_in_3way = ar_in_3way[ar_in_3way > 0]
    
    etc_in_3way = np.unique(etc_mask.values[threeway_mask > 0])
    etc_in_3way = etc_in_3way[etc_in_3way > 0]
    
    validated_triplets = []
    triplet_details = {}
    valid_mcs = []
    valid_ar = []
    valid_etc = []
    
    # For each potential triplet, verify mutual overlap and thresholds
    for mcs_id in mcs_in_3way:
        for ar_id in ar_in_3way:
            for etc_id in etc_in_3way:
                # Get pixels for each track
                mcs_pixels = (mcs_mask == mcs_id)
                ar_pixels = (ar_mask == ar_id)
                etc_pixels = (etc_mask == etc_id)
                
                # Find the intersection of all three tracks in 3-way regions
                mutual_overlap = mcs_pixels & ar_pixels & etc_pixels & (threeway_mask > 0)
                overlap_count = mutual_overlap.sum().compute().item()
                
                if overlap_count > 0:
                    # Calculate overlap fractions for each track
                    mcs_total = mcs_pixels.sum().compute().item()
                    ar_total = ar_pixels.sum().compute().item()
                    etc_total = etc_pixels.sum().compute().item()
                    
                    mcs_fraction = overlap_count / mcs_total if mcs_total > 0 else 0
                    ar_fraction = overlap_count / ar_total if ar_total > 0 else 0
                    etc_fraction = overlap_count / etc_total if etc_total > 0 else 0
                    
                    # Check if all tracks meet their thresholds
                    if (mcs_fraction >= mcs_thresh and 
                        ar_fraction >= ar_thresh and 
                        etc_fraction >= etc_thresh):
                        
                        triplet = (mcs_id, ar_id, etc_id)
                        validated_triplets.append(triplet)
                        valid_mcs.append(mcs_id)
                        valid_ar.append(ar_id)
                        valid_etc.append(etc_id)
                        
                        triplet_details[triplet] = {
                            'overlap_pixels': overlap_count,
                            'mcs_fraction': mcs_fraction,
                            'ar_fraction': ar_fraction,
                            'etc_fraction': etc_fraction
                        }
    
    # Remove duplicates
    valid_mcs = np.unique(valid_mcs)
    valid_ar = np.unique(valid_ar)
    valid_etc = np.unique(valid_etc)
    
    if verbose:
        print(f"Validated 3-way overlaps:")
        print(f"  {len(validated_triplets)} triplets found")
        print(f"  MCS tracks: {len(valid_mcs)} - {valid_mcs}")
        print(f"  AR tracks: {len(valid_ar)} - {valid_ar}")
        print(f"  ETC tracks: {len(valid_etc)} - {valid_etc}")
        if len(validated_triplets) > 0 and len(validated_triplets) <= 5:
            print(f"  Triplets: {validated_triplets}")
    
    return {
        'mcs_tracks': valid_mcs,
        'ar_tracks': valid_ar,
        'etc_tracks': valid_etc,
        'validated_triplets': validated_triplets,
        'triplet_details': triplet_details
    }


def filter_mcs_tc_overlaps(mcs_mask, tc_mask, overlap_threshold=0.15, verbose=True):
    """
    Filter MCS tracks that significantly overlap with tropical cyclones using vectorized operations.
    
    This function follows the efficient approach from find_overlapping_tracks_and_pairs:
    - Creates binary masks and sums them to find overlap regions
    - Extracts MCS track IDs from overlap regions using vectorized operations
    - Calculates overlap fractions efficiently without loops
    - Returns filtered MCS mask with TC-overlapping tracks removed
    
    Parameters:
    -----------
    mcs_mask : xarray.DataArray
        Original MCS mask with track IDs (values > 0 indicate tracks)
    tc_mask : xarray.DataArray  
        TC mask with track IDs (values > 0 indicate tracks)
    overlap_threshold : float, default=0.15
        Fractional overlap threshold (0.15 = 15%). MCS tracks with overlap 
        fraction >= threshold will be identified for removal
    verbose : bool, default=True
        Whether to print progress information
        
    Returns:
    --------
    dict : Results containing:
        - 'mcs_filtered': Filtered MCS mask with TC-overlapping tracks removed
        - 'removed_tracks': Boolean mask of pixels belonging to removed tracks
        - 'summary': Summary statistics
    """

    # Create binary masks (1 where feature exists, 0 elsewhere)
    mcs_binary = xr.where(mcs_mask > 0, 1, 0)
    tc_binary = xr.where(tc_mask > 0, 1, 0)
    # Sum binary masks to identify overlap regions
    binary_sum_mask = mcs_binary + tc_binary

    # Find overlap condition: pixels where both features exist AND no NaN values
    overlap_condition = (
        (binary_sum_mask >= 2) &  # Both features present
        ~np.isnan(binary_sum_mask) &  # No NaN in sum
        ~np.isnan(mcs_mask) &  # No NaN in MCS mask
        ~np.isnan(tc_mask)  # No NaN in TC mask
    )

    # Check if any overlaps exist
    total_overlap_pixels = overlap_condition.sum().compute().item()
    if total_overlap_pixels == 0:
        if verbose:
            print("No MCS-TC overlaps found, returning original mask")
        return {
            'mcs_filtered': mcs_mask,
            'removed_tracks': xr.zeros_like(mcs_mask, dtype=bool),
            'summary': {
                'initial_count': len(np.unique(mcs_mask.values[mcs_mask.values > 0])),
                'removed_count': 0,
                'final_count': len(np.unique(mcs_mask.values[mcs_mask.values > 0])),
                'retention_rate': 1.0
            }
        }
    
    if verbose:
        print(f"Found {total_overlap_pixels} overlap pixels")

    # Extract MCS track IDs that appear in overlap regions (vectorized operation)
    tracks_in_overlap, overlap_counts = np.unique(
        mcs_mask.where(overlap_condition), return_counts=True
    )

    # Remove NaN values from the results
    valid_mask = ~np.isnan(tracks_in_overlap)
    tracks_in_overlap = tracks_in_overlap[valid_mask]
    overlap_counts = overlap_counts[valid_mask]
    
    # Edge case: No valid tracks found in overlap regions
    if len(tracks_in_overlap) == 0:
        if verbose:
            print("No valid MCS tracks found in overlap regions, returning original mask")
        return {
            'mcs_filtered': mcs_mask,
            'removed_tracks': xr.zeros_like(mcs_mask, dtype=bool),
            'summary': {
                'initial_count': len(np.unique(mcs_mask.values[mcs_mask.values > 0])),
                'removed_count': 0,
                'final_count': len(np.unique(mcs_mask.values[mcs_mask.values > 0])),
                'retention_rate': 1.0
            }
        }
    
    if verbose:
        print(f"Found {len(tracks_in_overlap)} MCS tracks with TC overlap")

    # Get total pixel counts for ALL MCS tracks (not just overlapping ones)
    unique_tracks, total_counts = np.unique(mcs_mask, return_counts=True)
    
    # Create lookup dictionary for efficient access
    track_counts_dict = {}
    for track_id, count in zip(unique_tracks, total_counts):
        if not np.isnan(track_id):
            track_counts_dict[track_id] = count

    # Get total pixel counts for tracks that have overlaps (vectorized lookup)
    total_pixel_counts = np.array([
        track_counts_dict.get(track, 0) for track in tracks_in_overlap
    ])

    # Check for zero-sized tracks (should not happen but safety check)
    if np.any(total_pixel_counts == 0):
        if verbose:
            print("Warning: Found tracks with zero pixels, filtering them out")
        valid_indices = total_pixel_counts > 0
        tracks_in_overlap = tracks_in_overlap[valid_indices]
        overlap_counts = overlap_counts[valid_indices]
        total_pixel_counts = total_pixel_counts[valid_indices]
        
        if len(tracks_in_overlap) == 0:
            if verbose:
                print("No valid tracks remaining after filtering, returning original mask")
            return {
                'mcs_filtered': mcs_mask,
                'removed_tracks': xr.zeros_like(mcs_mask, dtype=bool),
                'summary': {
                    'initial_count': len(np.unique(mcs_mask.values[mcs_mask.values > 0])),
                    'removed_count': 0,
                    'final_count': len(np.unique(mcs_mask.values[mcs_mask.values > 0])),
                    'retention_rate': 1.0
                }
            }
    
    # Calculate overlap fractions for all tracks at once (vectorized)
    overlap_fractions = overlap_counts / total_pixel_counts

    # Find tracks that exceed the overlap threshold (vectorized comparison)
    threshold_mask = overlap_fractions >= overlap_threshold
    tracks_exceeding = tracks_in_overlap[threshold_mask]
    
    if verbose:
        n_exceeding = len(tracks_exceeding)
        n_total = len(tracks_in_overlap)
        print(f"Tracks exceeding {overlap_threshold:.1%} threshold: {n_exceeding}/{n_total}")
        
        if n_exceeding > 0 and n_exceeding <= 10:  # Show details for small numbers
            exceeding_fractions = overlap_fractions[threshold_mask]
            for track_id, fraction in zip(tracks_exceeding, exceeding_fractions):
                print(f"  Track {int(track_id)}: {fraction:.0%} overlap")
    
    # Edge case: No tracks exceed threshold
    if len(tracks_exceeding) == 0:
        if verbose:
            print("No tracks exceed overlap threshold, returning original mask")
        return {
            'mcs_filtered': mcs_mask,
            'removed_tracks': xr.zeros_like(mcs_mask, dtype=bool),
            'summary': {
                'initial_count': len(np.unique(mcs_mask.values[mcs_mask.values > 0])),
                'removed_count': 0,
                'final_count': len(np.unique(mcs_mask.values[mcs_mask.values > 0])),
                'retention_rate': 1.0
            }
        }

    # Create filtered mask using vectorized operations
    mcs_filtered = mcs_mask.copy()
    
    # Remove tracks that exceed the overlap threshold (vectorized operation)
    tracks_to_remove_mask = np.isin(mcs_mask, tracks_exceeding)
    
    # Set pixels belonging to tracks that exceed threshold to 0
    mcs_filtered = xr.where(tracks_to_remove_mask, 0, mcs_filtered)
    
    # Compile detailed overlap statistics
    original_tracks = np.unique(mcs_mask.values[mcs_mask.values > 0])
    remaining_tracks = np.unique(mcs_filtered.values[mcs_filtered.values > 0])

    n_mcs_original = len(original_tracks)
    n_mcs_final = len(remaining_tracks)
    n_removed = n_mcs_original - n_mcs_final
    retention_rate = n_mcs_final / n_mcs_original if n_mcs_original > 0 else 0.0
    
    if verbose:
        print(f"\nTC Filtering Results:")
        print(f"  Original MCS tracks: {n_mcs_original}")
        print(f"  MCS tracks removed due to TC overlap: {n_removed}")
        print(f"  Remaining tracks: {n_mcs_final}")
        print(f"  Retention rate: {retention_rate:.1%}")
    
    # Prepare summary
    summary = {
        'initial_count': n_mcs_original,
        'removed_count': n_removed,
        'final_count': n_mcs_final,
        'retention_rate': retention_rate
    }
    
    return {
        'mcs_filtered': mcs_filtered,
        'removed_tracks': tracks_to_remove_mask,
        'summary': summary
    }


def promote_dual_etc_overlaps_to_3way(mcs_ar_pairs_2way, ar_etc_pairs_2way, mcs_etc_pairs_2way, 
                                     existing_3way_tracks=None, verbose=True):
    """
    Identifies ETC tracks with dual 2-way overlaps and promotes them to 3-way,
    including ALL tracks that have 2-way overlaps with any of the promoted tracks.
    
    This ensures comprehensive grouping - if AR 184 gets promoted to 3-way with MCS 16470 and ETC 142,
    then ALL MCS tracks that overlap with AR 184 (e.g., MCS 16443) are also included in 3-way.
    
    Parameters:
    -----------
    mcs_ar_pairs_2way : list
        List of (mcs_track, ar_track) tuples for MCS-AR 2-way pairs
    ar_etc_pairs_2way : list
        List of (ar_track, etc_track) tuples for AR-ETC 2-way pairs  
    mcs_etc_pairs_2way : list
        List of (mcs_track, etc_track) tuples for MCS-ETC 2-way pairs
    existing_3way_tracks : dict, optional
        Existing 3-way tracks to merge with new ones
    verbose : bool, default=True
        Whether to print detailed information
        
    Returns:
    --------
    dict : Results with comprehensive track grouping
    """
    
    if verbose:
        print("🔍 IDENTIFYING ETC TRACKS WITH DUAL 2-WAY OVERLAPS")
        print("="*55)
    
    # Step 1: Find ETC tracks with dual overlaps (same as before)
    etc_tracks_in_ar_etc = set([pair[1] for pair in ar_etc_pairs_2way])
    etc_tracks_in_mcs_etc = set([pair[1] for pair in mcs_etc_pairs_2way])
    dual_etc_tracks = etc_tracks_in_ar_etc.intersection(etc_tracks_in_mcs_etc)
    
    if verbose:
        print(f"ETC tracks in AR-ETC pairs: {sorted(etc_tracks_in_ar_etc)}")
        print(f"ETC tracks in MCS-ETC pairs: {sorted(etc_tracks_in_mcs_etc)}")
        print(f"ETC tracks with dual overlaps: {sorted(dual_etc_tracks)}")
    
    if not dual_etc_tracks:
        if verbose:
            print("No ETC tracks found with dual 2-way overlaps.")
        
        return {
            'promoted_triplets': [],
            'updated_3way_tracks': existing_3way_tracks or {'mcs_tracks': np.array([]), 'ar_tracks': np.array([]), 'etc_tracks': np.array([])},
            'updated_2way_pairs': {
                'mcs_ar_pairs': mcs_ar_pairs_2way,
                'ar_etc_pairs': ar_etc_pairs_2way,
                'mcs_etc_pairs': mcs_etc_pairs_2way
            },
            'promotion_details': {}
        }
    
    # Step 2: Create initial triplets from dual ETC overlaps
    promoted_triplets = []
    core_promoted_tracks = {'mcs': set(), 'ar': set(), 'etc': set()}
    
    for etc_track in dual_etc_tracks:
        if verbose:
            print(f"\n🎯 Processing ETC track {etc_track}:")
        
        mcs_partners = [pair[0] for pair in mcs_etc_pairs_2way if pair[1] == etc_track]
        ar_partners = [pair[0] for pair in ar_etc_pairs_2way if pair[1] == etc_track]
        
        if verbose:
            print(f"  MCS partners: {mcs_partners}")
            print(f"  AR partners: {ar_partners}")
        
        for mcs_track in mcs_partners:
            for ar_track in ar_partners:
                triplet = (mcs_track, ar_track, etc_track)
                promoted_triplets.append(triplet)
                core_promoted_tracks['mcs'].add(mcs_track)
                core_promoted_tracks['ar'].add(ar_track)
                core_promoted_tracks['etc'].add(etc_track)
                
                if verbose:
                    print(f"  ✅ Core triplet: MCS {mcs_track} - AR {ar_track} - ETC {etc_track}")
    
    # Step 3: COMPREHENSIVE EXPANSION - Find ALL tracks with 2-way overlaps to promoted tracks
    comprehensive_3way_tracks = {'mcs': set(), 'ar': set(), 'etc': set()}
    
    # Start with core promoted tracks
    comprehensive_3way_tracks['mcs'].update(core_promoted_tracks['mcs'])
    comprehensive_3way_tracks['ar'].update(core_promoted_tracks['ar'])
    comprehensive_3way_tracks['etc'].update(core_promoted_tracks['etc'])
    
    if verbose:
        print(f"\n🌟 COMPREHENSIVE EXPANSION:")
        print(f"Core promoted tracks: MCS={sorted(core_promoted_tracks['mcs'])}, AR={sorted(core_promoted_tracks['ar'])}, ETC={sorted(core_promoted_tracks['etc'])}")
    
    # Add ALL MCS tracks that overlap with promoted AR tracks
    for mcs_track, ar_track in mcs_ar_pairs_2way:
        if ar_track in core_promoted_tracks['ar']:
            comprehensive_3way_tracks['mcs'].add(mcs_track)
            if verbose and mcs_track not in core_promoted_tracks['mcs']:
                print(f"  + Added MCS {mcs_track} (overlaps with promoted AR {ar_track})")
    
    # Add ALL AR tracks that overlap with promoted MCS tracks
    for mcs_track, ar_track in mcs_ar_pairs_2way:
        if mcs_track in core_promoted_tracks['mcs']:
            comprehensive_3way_tracks['ar'].add(ar_track)
            if verbose and ar_track not in core_promoted_tracks['ar']:
                print(f"  + Added AR {ar_track} (overlaps with promoted MCS {mcs_track})")
    
    # Add ALL AR tracks that overlap with promoted ETC tracks
    for ar_track, etc_track in ar_etc_pairs_2way:
        if etc_track in core_promoted_tracks['etc']:
            comprehensive_3way_tracks['ar'].add(ar_track)
            if verbose and ar_track not in core_promoted_tracks['ar']:
                print(f"  + Added AR {ar_track} (overlaps with promoted ETC {etc_track})")
    
    # Add ALL ETC tracks that overlap with promoted AR tracks
    for ar_track, etc_track in ar_etc_pairs_2way:
        if ar_track in core_promoted_tracks['ar']:
            comprehensive_3way_tracks['etc'].add(etc_track)
            if verbose and etc_track not in core_promoted_tracks['etc']:
                print(f"  + Added ETC {etc_track} (overlaps with promoted AR {ar_track})")
    
    # Add ALL MCS tracks that overlap with promoted ETC tracks
    for mcs_track, etc_track in mcs_etc_pairs_2way:
        if etc_track in core_promoted_tracks['etc']:
            comprehensive_3way_tracks['mcs'].add(mcs_track)
            if verbose and mcs_track not in core_promoted_tracks['mcs']:
                print(f"  + Added MCS {mcs_track} (overlaps with promoted ETC {etc_track})")
    
    # Add ALL ETC tracks that overlap with promoted MCS tracks
    for mcs_track, etc_track in mcs_etc_pairs_2way:
        if mcs_track in core_promoted_tracks['mcs']:
            comprehensive_3way_tracks['etc'].add(etc_track)
            if verbose and etc_track not in core_promoted_tracks['etc']:
                print(f"  + Added ETC {etc_track} (overlaps with promoted MCS {mcs_track})")
    
    # Step 4: Remove ALL comprehensive tracks from 2-way pairs
    updated_mcs_ar_pairs = [
        pair for pair in mcs_ar_pairs_2way 
        if pair[0] not in comprehensive_3way_tracks['mcs'] and pair[1] not in comprehensive_3way_tracks['ar']
    ]
    
    updated_ar_etc_pairs = [
        pair for pair in ar_etc_pairs_2way 
        if pair[0] not in comprehensive_3way_tracks['ar'] and pair[1] not in comprehensive_3way_tracks['etc']
    ]
    
    updated_mcs_etc_pairs = [
        pair for pair in mcs_etc_pairs_2way 
        if pair[0] not in comprehensive_3way_tracks['mcs'] and pair[1] not in comprehensive_3way_tracks['etc']
    ]
    
    # Convert to arrays and merge with existing 3-way tracks
    comprehensive_mcs = np.array(sorted(comprehensive_3way_tracks['mcs']))
    comprehensive_ar = np.array(sorted(comprehensive_3way_tracks['ar']))
    comprehensive_etc = np.array(sorted(comprehensive_3way_tracks['etc']))
    
    if existing_3way_tracks:
        final_mcs_3way = np.unique(np.concatenate([existing_3way_tracks['mcs_tracks'], comprehensive_mcs]))
        final_ar_3way = np.unique(np.concatenate([existing_3way_tracks['ar_tracks'], comprehensive_ar]))
        final_etc_3way = np.unique(np.concatenate([existing_3way_tracks['etc_tracks'], comprehensive_etc]))
    else:
        final_mcs_3way = comprehensive_mcs
        final_ar_3way = comprehensive_ar
        final_etc_3way = comprehensive_etc
    
    if verbose:
        print(f"\n📊 COMPREHENSIVE PROMOTION SUMMARY:")
        print(f"  Core triplets promoted: {len(promoted_triplets)}")
        print(f"  Comprehensive 3-way tracks:")
        print(f"    MCS: {len(final_mcs_3way)} tracks {sorted(final_mcs_3way)}")
        print(f"    AR: {len(final_ar_3way)} tracks {sorted(final_ar_3way)}")
        print(f"    ETC: {len(final_etc_3way)} tracks {sorted(final_etc_3way)}")
        print(f"  Updated 2-way pairs:")
        print(f"    MCS-AR: {len(updated_mcs_ar_pairs)} (was {len(mcs_ar_pairs_2way)})")
        print(f"    AR-ETC: {len(updated_ar_etc_pairs)} (was {len(ar_etc_pairs_2way)})")
        print(f"    MCS-ETC: {len(updated_mcs_etc_pairs)} (was {len(mcs_etc_pairs_2way)})")
    
    return {
        'promoted_triplets': promoted_triplets,
        'updated_3way_tracks': {
            'mcs_tracks': final_mcs_3way,
            'ar_tracks': final_ar_3way,
            'etc_tracks': final_etc_3way
        },
        'updated_2way_pairs': {
            'mcs_ar_pairs': updated_mcs_ar_pairs,
            'ar_etc_pairs': updated_ar_etc_pairs,
            'mcs_etc_pairs': updated_mcs_etc_pairs
        },
        'promotion_details': {
            'core_promoted_tracks': core_promoted_tracks,
            'comprehensive_tracks': comprehensive_3way_tracks
        }
    }


def process_single_timestep_overlaps(_ds, verbose=True):
    """
    Process co-occurrence feature overlaps for a single time step.
    
    Parameters:
    -----------
    _ds : xarray.Dataset
        Single time step dataset containing mcs_mask, ar_mask, etc_mask, tc_mask
    verbose : bool, default=True
        Whether to print progress information
        
    Returns:
    --------
    dict : Dictionary containing all isolated and overlap masks for this time step
    """
    
    if verbose:
        print(f"Processing single time step overlap analysis...")
    
    # Define thresholds based on atmospheric scale hierarchy
    mcs_thresh = 0.20  # 20% - smallest features, need substantial overlap
    ar_thresh = 0.10   # 10% - medium features
    etc_thresh_2way = 0.05  # 5% - largest features for 2-way overlaps
    etc_thresh_3way = 0.00  # 0% - environmental context for 3-way overlaps
    
    # ===== STEP 1: MCS-TC FILTERING =====
    if verbose:
        print("  Step 1: Filtering MCS-TC overlaps...")
    
    tc_filtering_results = filter_mcs_tc_overlaps(
        mcs_mask=_ds.mcs_mask,
        tc_mask=_ds.tc_mask,
        overlap_threshold=0.10,  # 10% threshold
        verbose=verbose
    )
    
    mcs_filtered = tc_filtering_results['mcs_filtered']
    
    # ===== STEP 2: CREATE BINARY MASKS =====
    if verbose:
        print("  Step 2: Creating binary masks...")
        
    # Create binary masks (1 where feature exists, 0 elsewhere)
    mcs_binary = xr.where(mcs_filtered > 0, 1, 0)
    ar_binary = xr.where(_ds.ar_mask > 0, 1, 0)
    etc_binary = xr.where(_ds.etc_mask > 0, 1, 0)

    # Create pairwise summation masks
    mcs_ar_sum = mcs_binary + ar_binary
    ar_etc_sum = ar_binary + etc_binary  
    mcs_etc_sum = mcs_binary + etc_binary

    # Create 3-way summation mask
    mcs_ar_etc_sum = mcs_binary + ar_binary + etc_binary
    
    # ===== STEP 3: 3-WAY OVERLAPS =====
    if verbose:
        print("  Step 3: Finding 3-way overlaps...")
    
    validated_3way = find_true_3way_overlaps(
        mcs_filtered, _ds.ar_mask, _ds.etc_mask, mcs_ar_etc_sum,
        mcs_thresh, ar_thresh, etc_thresh_3way, verbose=verbose
    )

    # Extract validated tracks
    mcs_tracks_with_ar_etc_overlap = validated_3way['mcs_tracks']
    ar_tracks_with_mcs_etc_overlap = validated_3way['ar_tracks']
    etc_tracks_with_mcs_ar_overlap = validated_3way['etc_tracks']
    validated_triplets = validated_3way['validated_triplets']
    
    # ===== STEP 4: 2-WAY PAIRS =====
    if verbose:
        print("  Step 4: Finding 2-way pairs...")
    
    # Find all potential 2-way overlaps
    mcs_ar_pairs_all = find_overlapping_tracks_and_pairs(mcs_filtered, _ds.ar_mask, mcs_ar_sum,
                                                        mcs_thresh, 0.0, overlap_threshold=1,
                                                        feature1_name="MCS", feature2_name="AR", verbose=verbose)

    ar_etc_pairs_all = find_overlapping_tracks_and_pairs(_ds.ar_mask, _ds.etc_mask, ar_etc_sum,
                                                        ar_thresh, 0.01, overlap_threshold=1,
                                                        feature1_name="AR", feature2_name="ETC", verbose=verbose)

    mcs_etc_pairs_all = find_overlapping_tracks_and_pairs(mcs_filtered, _ds.etc_mask, mcs_etc_sum,
                                                         mcs_thresh, 0.0, overlap_threshold=1,
                                                         feature1_name="MCS", feature2_name="ETC", verbose=verbose)

    # Filter out tracks that ACTUALLY EXCEED 3-way thresholds (not just participate)
    mcs_ar_pairs_2way_only = [pair for pair in mcs_ar_pairs_all['valid_pairs']
                             if pair[0] not in mcs_tracks_with_ar_etc_overlap 
                             and pair[1] not in ar_tracks_with_mcs_etc_overlap]

    ar_etc_pairs_2way_only = [pair for pair in ar_etc_pairs_all['valid_pairs']
                             if pair[0] not in ar_tracks_with_mcs_etc_overlap 
                             and pair[1] not in etc_tracks_with_mcs_ar_overlap]

    mcs_etc_pairs_2way_only = [pair for pair in mcs_etc_pairs_all['valid_pairs']
                              if pair[0] not in mcs_tracks_with_ar_etc_overlap 
                              and pair[1] not in etc_tracks_with_mcs_ar_overlap]
    
    # ===== STEP 5: DUAL ETC OVERLAP PROMOTION =====
    if verbose:
        print("  Step 5: Promoting dual ETC overlaps...")
        
    comprehensive_results = promote_dual_etc_overlaps_to_3way(
        mcs_ar_pairs_2way_only, ar_etc_pairs_2way_only, mcs_etc_pairs_2way_only,
        existing_3way_tracks={
            'mcs_tracks': mcs_tracks_with_ar_etc_overlap,
            'ar_tracks': ar_tracks_with_mcs_etc_overlap,
            'etc_tracks': etc_tracks_with_mcs_ar_overlap
        },
        verbose=verbose
    )

    # Update with comprehensive results
    mcs_ar_pairs_2way_only = comprehensive_results['updated_2way_pairs']['mcs_ar_pairs']
    ar_etc_pairs_2way_only = comprehensive_results['updated_2way_pairs']['ar_etc_pairs']
    mcs_etc_pairs_2way_only = comprehensive_results['updated_2way_pairs']['mcs_etc_pairs']

    # Update 3-way tracks with comprehensive grouping
    mcs_tracks_with_ar_etc_overlap = comprehensive_results['updated_3way_tracks']['mcs_tracks']
    ar_tracks_with_mcs_etc_overlap = comprehensive_results['updated_3way_tracks']['ar_tracks']
    etc_tracks_with_mcs_ar_overlap = comprehensive_results['updated_3way_tracks']['etc_tracks']
    
    # ===== STEP 6: CREATE VISUALIZATION MASKS =====
    if verbose:
        print("  Step 6: Creating output masks...")
    
    # Extract track lists from 2-way pairs
    if mcs_ar_pairs_2way_only:
        mcs_tracks_paired_with_ar = np.unique([pair[0] for pair in mcs_ar_pairs_2way_only])
        ar_tracks_paired_with_mcs = np.unique([pair[1] for pair in mcs_ar_pairs_2way_only])
    else:
        mcs_tracks_paired_with_ar = np.array([])
        ar_tracks_paired_with_mcs = np.array([])

    if ar_etc_pairs_2way_only:
        ar_tracks_paired_with_etc = np.unique([pair[0] for pair in ar_etc_pairs_2way_only])
        etc_tracks_paired_with_ar = np.unique([pair[1] for pair in ar_etc_pairs_2way_only])
    else:
        ar_tracks_paired_with_etc = np.array([])
        etc_tracks_paired_with_ar = np.array([])

    if mcs_etc_pairs_2way_only:
        mcs_tracks_paired_with_etc = np.unique([pair[0] for pair in mcs_etc_pairs_2way_only])
        etc_tracks_paired_with_mcs = np.unique([pair[1] for pair in mcs_etc_pairs_2way_only])
    else:
        mcs_tracks_paired_with_etc = np.array([])
        etc_tracks_paired_with_mcs = np.array([])

    # Create 2-way pair masks
    mcs_ar_overlap_mask_pairs = xr.where(
        mcs_filtered.isin(mcs_tracks_paired_with_ar), 
        mcs_filtered, 
        0
    )
    ar_mcs_overlap_mask_pairs = xr.where(
        _ds.ar_mask.isin(ar_tracks_paired_with_mcs), 
        _ds.ar_mask, 
        0
    )

    ar_etc_overlap_mask_pairs = xr.where(
        _ds.ar_mask.isin(ar_tracks_paired_with_etc), 
        _ds.ar_mask, 
        0
    )
    etc_ar_overlap_mask_pairs = xr.where(
        _ds.etc_mask.isin(etc_tracks_paired_with_ar), 
        _ds.etc_mask, 
        0
    )

    mcs_etc_overlap_mask_pairs = xr.where(
        mcs_filtered.isin(mcs_tracks_paired_with_etc), 
        mcs_filtered, 
        0
    )
    etc_mcs_overlap_mask_pairs = xr.where(
        _ds.etc_mask.isin(etc_tracks_paired_with_mcs), 
        _ds.etc_mask, 
        0
    )

    # Create 3-way overlap masks
    mcs_ar_etc_overlap_mask = xr.where(
        mcs_filtered.isin(mcs_tracks_with_ar_etc_overlap), 
        mcs_filtered, 
        0
    )
    ar_mcs_etc_overlap_mask = xr.where(
        _ds.ar_mask.isin(ar_tracks_with_mcs_etc_overlap), 
        _ds.ar_mask, 
        0
    )
    etc_mcs_ar_overlap_mask = xr.where(
        _ds.etc_mask.isin(etc_tracks_with_mcs_ar_overlap), 
        _ds.etc_mask, 
        0
    )

    # Combine all MCS, AR, ETC tracks involved in any overlaps
    mcs_tracks_all_overlap = np.concatenate((mcs_tracks_paired_with_ar, mcs_tracks_paired_with_etc, mcs_tracks_with_ar_etc_overlap))
    ar_tracks_all_overlap = np.concatenate((ar_tracks_paired_with_mcs, ar_tracks_paired_with_etc, ar_tracks_with_mcs_etc_overlap))
    etc_tracks_all_overlap = np.concatenate((etc_tracks_paired_with_mcs, etc_tracks_paired_with_ar, etc_tracks_with_mcs_ar_overlap))

    # Create isolated masks excluding all overlapping tracks
    mcs_isolated_mask = xr.where(
        mcs_filtered.isin(mcs_tracks_all_overlap), 
        0, 
        mcs_filtered,
    ).fillna(0)

    ar_isolated_mask = xr.where(
        _ds.ar_mask.isin(ar_tracks_all_overlap), 
        0,
        _ds.ar_mask,
    ).fillna(0)

    etc_isolated_mask = xr.where(
        _ds.etc_mask.isin(etc_tracks_all_overlap), 
        0, 
        _ds.etc_mask,
    ).fillna(0)

    if verbose:
        print(f"  ✅ Completed processing for this time step")
        print(f"    2-way pairs: MCS-AR={len(mcs_ar_pairs_2way_only)}, AR-ETC={len(ar_etc_pairs_2way_only)}, MCS-ETC={len(mcs_etc_pairs_2way_only)}")
        print(f"    3-way tracks: MCS={len(mcs_tracks_with_ar_etc_overlap)}, AR={len(ar_tracks_with_mcs_etc_overlap)}, ETC={len(etc_tracks_with_mcs_ar_overlap)}")

    # Return all masks and metadata
    return {
        # Original masks (TC-filtered for MCS)
        'mcs_mask': mcs_filtered,
        'ar_mask': _ds.ar_mask,
        'etc_mask': _ds.etc_mask,
        'tc_mask': _ds.tc_mask,
        'ccs_mask': _ds.ccs_mask,
        
        # Isolated masks
        'mcs_isolated_mask': mcs_isolated_mask,
        'ar_isolated_mask': ar_isolated_mask,
        'etc_isolated_mask': etc_isolated_mask,
        
        # 2-way overlap masks
        'mcs_ar_overlap_mask': mcs_ar_overlap_mask_pairs,
        'ar_mcs_overlap_mask': ar_mcs_overlap_mask_pairs,
        'ar_etc_overlap_mask': ar_etc_overlap_mask_pairs,
        'etc_ar_overlap_mask': etc_ar_overlap_mask_pairs,
        'mcs_etc_overlap_mask': mcs_etc_overlap_mask_pairs,
        'etc_mcs_overlap_mask': etc_mcs_overlap_mask_pairs,
        
        # 3-way overlap masks
        'mcs_ar_etc_overlap_mask': mcs_ar_etc_overlap_mask,
        'ar_mcs_etc_overlap_mask': ar_mcs_etc_overlap_mask,
        'etc_mcs_ar_overlap_mask': etc_mcs_ar_overlap_mask,
        
        # Metadata
        'mcs_ar_pairs_2way': mcs_ar_pairs_2way_only,
        'ar_etc_pairs_2way': ar_etc_pairs_2way_only,
        'mcs_etc_pairs_2way': mcs_etc_pairs_2way_only,
        'validated_triplets': validated_triplets,
        'tc_filtering_summary': tc_filtering_results['summary'],
        
        # Track lists for reference
        'mcs_tracks_3way': mcs_tracks_with_ar_etc_overlap,
        'ar_tracks_3way': ar_tracks_with_mcs_etc_overlap,
        'etc_tracks_3way': etc_tracks_with_mcs_ar_overlap
    }


def main():
    """
    Main function that processes the full time series dataset.
    """
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process co-occurrence feature overlaps')
    parser.add_argument('--source', type=str, required=True, help='Source name')
    parser.add_argument('--parallel', action='store_true', default=True,
                       help='Use parallel processing with Dask (default: True)')
    parser.add_argument('--no-parallel', action='store_false', dest='parallel',
                       help='Disable parallel processing')
    parser.add_argument('--workers', type=int, default=32,
                       help='Number of Dask workers (default: 32)')
    parser.add_argument('--threads-per-worker', type=int, default=1,
                       help='Number of threads per worker (default: 1)')
    parser.add_argument('--test-steps', type=int, default=None,
                       help='Number of time steps to process for testing (default: all)')
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    # Configuration from arguments
    source_name = args.source
    root_dir = "/pscratch/sd/w/wcmca1/hackathon/all_masks/"
    in_dir = f"{root_dir}/{source_name}_allmasks_hp8_v1.zarr"
    output_dir = "/pscratch/sd/w/wcmca1/hackathon/cof_masks/"
    output_path = f"{output_dir}/{source_name}_cofmasks_hp8_v1.zarr"
    
    # Parallel processing configuration
    parallel = args.parallel
    n_workers = args.workers
    threads_per_worker = args.threads_per_worker
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*80)
    print("CO-OCCURRENCE FEATURE OVERLAP PROCESSING")
    print("="*80)
    print(f"Source: {source_name}")
    print(f"Input: {in_dir}")
    print(f"Output: {output_dir}")
    print(f"Parallel processing: {parallel}")
    if parallel:
        print(f"Workers: {n_workers}, Threads per worker: {threads_per_worker}")
    
    # Setup Dask client
    client = setup_dask_client(parallel=parallel, n_workers=n_workers, threads_per_worker=threads_per_worker, logger=logger)
    
    try:
        # Load the full dataset
        print(f"\nLoading full dataset...")
        try:
            ds = xr.open_dataset(in_dir)
            # ds = ds.pipe(egh.attach_coords)  # Commented out for testing
            print(f"  ✅ Dataset loaded successfully")
            print(f"  Time steps: {len(ds.time)}")
            print(f"  Data variables: {list(ds.data_vars)}")
            print(f"  Spatial dimensions: {dict(ds.dims)}")
        except Exception as e:
            print(f"  ❌ Error loading dataset: {e}")
            return
        
        # Limit time steps for testing if requested
        if args.test_steps is not None:
            ds = ds.isel(time=slice(0, args.test_steps))
            print(f"  📋 Limited to {args.test_steps} time steps for testing")
        
        # Initialize streaming processing configuration
        print(f"\nInitializing streaming zarr processing...")
        time_coords = ds.time.values
        
        # Define all output variables
        mask_variables = [
            'mcs_mask', 'ar_mask', 'etc_mask', 'tc_mask', 'ccs_mask',
            'mcs_isolated_mask', 'ar_isolated_mask', 'etc_isolated_mask',
            'mcs_ar_overlap_mask', 'ar_mcs_overlap_mask',
            'ar_etc_overlap_mask', 'etc_ar_overlap_mask', 
            'mcs_etc_overlap_mask', 'etc_mcs_overlap_mask',
            'mcs_ar_etc_overlap_mask', 'ar_mcs_etc_overlap_mask', 'etc_mcs_ar_overlap_mask'
        ]
        
        # Add processing metadata
        attrs = ds.attrs.copy()
        attrs.update({
            'processing_info': 'Co-occurrence feature masks',
            'processing_date': str(np.datetime64('today')),
            'thresholds': 'MCS: 20%, AR: 10%, ETC: 5% (2-way), 0% (3-way)',
            'tc_filtering_threshold': '10%',
            'parallel_processing': str(parallel),
            'memory_approach': 'streaming'
        })
        if parallel:
            attrs['dask_workers'] = str(n_workers)
            attrs['threads_per_worker'] = str(threads_per_worker)
        
        # Stream processing and writing to zarr
        print(f"\nStreaming processing and writing {len(time_coords)} time steps to zarr...")
        
        # Use a simple default chunk size for time dimension
        # With zarr-path approach, serialization is minimal regardless of chunk size
        # Chunk size only affects processing efficiency and zarr I/O
        chunk_size_time = 48
        print(f"Using default chunk_size_time={chunk_size_time} for optimal processing and zarr I/O")
        
        try:
            # Initialize the zarr store structure (once only)
            logger.info("Initializing zarr store...")
            initialize_zarr_store(
                output_path=output_path,
                time_coords=time_coords,
                mask_variables=mask_variables,
                template_coords=ds.coords,
                attrs=attrs,
                chunk_size_time=chunk_size_time
            )
            
            # Stream process with chunked zarr writing
            successful_times = stream_process_to_zarr(
                ds=ds,
                time_coords=time_coords,
                mask_variables=mask_variables,
                output_path=output_path,
                template_coords=ds.coords,
                attrs=attrs,
                client=client,
                logger=logger,
                parallel=parallel,
                chunk_size_time=chunk_size_time,
                input_zarr_path=in_dir  # Pass the input zarr path for workers
            )
            
            logger.info(f"✅ Processing complete: {successful_times} time steps written to {output_path}")
            
        except Exception as e:
            logger.error(f"Error writing chunked zarr: {e}")
            print(f"  ❌ Error writing zarr: {e}")
            return
        
        # Store success info to print after Dask cleanup
        success_info = {
            'output_path': output_path,
            'successful_times': successful_times,
            'total_times': len(time_coords),
            'mask_variables': mask_variables,
            'parallel': parallel
        }

    finally:
        # Always cleanup client
        if client and parallel:
            # Suppress Dask shutdown messages by temporarily raising log level
            logging.getLogger('distributed').setLevel(logging.CRITICAL)
            logging.getLogger('distributed.worker').setLevel(logging.CRITICAL)
            logging.getLogger('distributed.nanny').setLevel(logging.CRITICAL)
            
            client.close()
        
        # Print success message after Dask cleanup (so it's always visible at the end)
        if 'success_info' in locals():
            print(f"\n{'='*80}")
            print(f"✅ PROCESSING COMPLETE!")
            print(f"{'='*80}")
            print(f"Output: {success_info['output_path']}")
            print(f"Time steps processed: {success_info['successful_times']}/{success_info['total_times']}")
            print(f"Variables created: {len(success_info['mask_variables'])}")
            print(f"Processing mode: {'Parallel' if success_info['parallel'] else 'Sequential'}")
            print(f"Memory approach: Streaming zarr writing")
            
            # Verify the saved dataset for sample statistics
            if success_info['successful_times'] > 0:
                try:
                    # Open the written zarr to get sample statistics
                    test_ds = xr.open_dataset(success_info['output_path'], engine='zarr')
                    sample_time = test_ds.time.values[0]
                    
                    print(f"\nSUMMARY STATISTICS:")
                    print(f"  Example from {sample_time}:")
                    print(f"    Dataset successfully written with {len(test_ds.time)} time steps")
                    print(f"    Variables created: {len([v for v in test_ds.data_vars if 'mask' in v])}")
                    print(f"    Spatial cells: {test_ds.sizes.get('cell', 'unknown')}")
                    test_ds.close()
                except Exception as e:
                    logger.warning(f"Could not read sample statistics: {e}")
                    print(f"  Dataset written but could not read sample statistics")
            print(f"{'='*80}\n")


if __name__ == "__main__":
    main()