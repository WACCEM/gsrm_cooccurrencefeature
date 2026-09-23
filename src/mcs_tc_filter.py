"""
Shared MCS-TC track-overlap filter.

This is the single definition of which MCS tracks are "TC-contaminated". It is used by two
pipeline steps, but only one of them excludes anything:

- scripts/make_mcs_swath_masks.py (Step 1): applies the exclusion to the native-cadence
  hourly MCS mask (overlap pooled over each aggregation window) *before* cloud-type
  classification (classify_cloud_types), so pixels freed by the exclusion get properly
  classified into DC/ND/ST/DZ from Tb/pr instead of being left at the zero cloud-type/
  cloud-precip values Step 1 assigns inside the MCS swath.
- scripts/make_cooccurrence_masks.py (Step 3): runs the same test on the aggregated 6-hourly
  swath for information only and removes nothing. The swath is a different mask from the
  hourly ones Step 1 tested (union footprint, coverage-priority pixel loss), so it flags
  borderline tracks (7-16% of 6-hourly frames in the September 2026 reprocessing) that Step 1
  correctly kept; removing them there would leave their precipitation in no category.

Historical context: this filter used to run only in Step 3, *after* Step 1 had already
zeroed cloud_types/dc_pr/st_pr/nd_pr/dz_pr for the full (pre-TC-filter) MCS swath, so any
track it removed had its precipitation attributed nowhere and fell through to "Residual" in
the downstream COF precipitation budget. Running it in Step 1, before the cloud-type
zeroing, closes that gap; keeping Step 3 informational avoids reopening it.

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import numpy as np
import xarray as xr

# MCS tracks with an overlap fraction against TC pixels at or above this threshold
# (within a single aggregation window) are excluded from the MCS mask entirely for that
# window by Step 1; Step 3 only uses it to count and log. Import this constant rather
# than hardcoding it.
MCS_TC_FILTER_THRESHOLD = 0.10


def filter_mcs_tc_overlaps(mcs_mask, tc_mask, overlap_threshold=MCS_TC_FILTER_THRESHOLD, verbose=True):
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
        MCS mask with track IDs (values > 0 indicate tracks). Any shape/dims are
        supported (e.g. a single aggregation window's (cell,) swath, as used by Step 3,
        or a full native-cadence (time, cell) window, as used by Step 1) - the overlap
        fraction is computed over all elements, flattened.
    tc_mask : xarray.DataArray
        TC mask with track IDs (values > 0 indicate tracks). Same shape as mcs_mask.
    overlap_threshold : float, default=MCS_TC_FILTER_THRESHOLD (0.10 = 10%)
        Fractional overlap threshold. MCS tracks with overlap fraction >= threshold
        will be identified for removal.
    verbose : bool, default=True
        Whether to print progress information.

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
