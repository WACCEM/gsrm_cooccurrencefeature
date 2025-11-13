#!/usr/bin/env python
"""
Plot all 4 feature masks (MCS, AR, ETC, TC) for a time range with Dask parallelization.

Usage:
python plot_feature_masks.py --source scream --start 2019-10-31 --end 2019-11-02 [options]

Optional arguments:
--parallel 0 (serial), 1 (parallel) - default: 1
--figdir output_directory - default: /global/cfs/cdirs/m1867/zfeng/hk25/quicklooks/
--figsize width height - figure size in inches, default: 12 6.75
--dpi DPI for output figures, default: 200
--workers number of dask workers, default: 4

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import argparse
import os
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib as mpl
# Set non-GUI backend for thread safety
mpl.use('agg')
import matplotlib.pyplot as plt
# Configure matplotlib to be more memory-efficient
plt.rcParams['figure.max_open_warning'] = 0  # Disable the warning
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import cartopy.crs as ccrs
import cartopy.feature as cf
import warnings
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from easygems import healpix as egh
from src.zarr_tools import setup_dask_client

# Suppress warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def plot_features_with_cooccurrence(ds, title="", figsize=(12, 6), 
                                   base_colors=None, base_alphas=None,
                                   cooccur_colors=None, cooccur_alphas=None,
                                   figname=None, dpi=200, 
                                   title_fontsize=14, legend_fontsize=12, 
                                   show_track_labels=True, track_fontsize=8):
    """
    Plot all 4 individual feature masks with co-occurring feature masks overlaid.
    Thread-safe version matching plot_all_feature_masks() styling.
    
    Co-occurring masks are created by combining overlap masks:
    - MCS-AR: union of mcs_ar_overlap_mask + ar_mcs_overlap_mask
    - MCS-ETC: union of mcs_etc_overlap_mask + etc_mcs_overlap_mask
    - AR-ETC: union of ar_etc_overlap_mask + etc_ar_overlap_mask
    - MCS-AR-ETC: union of mcs_ar_etc_overlap_mask + ar_mcs_etc_overlap_mask + etc_mcs_ar_overlap_mask
    
    Parameters:
    - ds: xarray Dataset containing the mask variables (assumes single time step)
    - title: optional title for the plot
    - figsize: figure size tuple
    - base_colors: dict with colors for individual features {'mcs', 'ar', 'etc', 'tc'}
    - base_alphas: dict with alpha values for individual features
    - cooccur_colors: dict with colors for co-occurring features {'mcs_ar', 'mcs_etc', 'ar_etc', 'mcs_ar_etc'}
    - cooccur_alphas: dict with alpha values for co-occurring features
    - figname: if provided, save figure as PNG with this filename
    - dpi: DPI for figure display and saving (default 200)
    - title_fontsize: font size for the plot title (default 14)
    - legend_fontsize: font size for the legend (default 12)
    - show_track_labels: whether to show track ID labels for co-occurring MCS, ETC and AR features (default True)
    - track_fontsize: font size for track labels (default 8)
    
    Returns:
    - figname: filename of the saved figure (or None if not saved)
    """
    # Set default colors for individual features (match plot_all_feature_masks)
    if base_colors is None:
        base_colors = {
            'mcs': 'lightskyblue',
            'ar': 'darkorange', 
            'etc': 'limegreen',
            'tc': 'crimson'
        }
    
    # Set default alphas for individual features (match plot_all_feature_masks)
    if base_alphas is None:
        # base_alphas = {
        #     'mcs': 0.8,
        #     'ar': 0.5,
        #     'etc': 0.3,
        #     'tc': 0.6
        # }
        # Reduced by half to make co-occurring features stand out more
        base_alphas = {
            'mcs': 0.4,
            'ar': 0.25,
            'etc': 0.15,
            'tc': 0.3
        }
    
    # Set default colors for co-occurring features (blend individual colors)
    if cooccur_colors is None:
        cooccur_colors = {
            'mcs_ar': '#FF7F50',        # Coral (blend of blue + orange)
            'mcs_etc': '#40E0D0',       # Turquoise (blend of blue + green)
            'ar_etc': '#9ACD32',        # Yellow-green (blend of orange + green)
            'mcs_ar_etc': '#C532C0'     # Dark purple (blend of all three)
        }
    
    # Set default alphas for co-occurring features
    if cooccur_alphas is None:
        cooccur_alphas = {
            'mcs_ar': 0.4,
            'mcs_etc': 0.4,
            'ar_etc': 0.4,
            'mcs_ar_etc': 0.4
        }
    
    # Create base map (match plot_all_feature_masks exactly)
    projection = ccrs.Robinson(central_longitude=-135)
    fig = plt.figure(figsize=figsize, dpi=dpi, facecolor='w')
    ax = plt.subplot(111, projection=projection)
    ax.set_global()
    ax.add_feature(cf.COASTLINE, linewidth=0.8)
    ax.add_feature(cf.BORDERS, linewidth=0.4)

    # Gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 10}
    gl.ylabel_style = {"size": 10}

    # Set title
    plot_title = title if title else "Feature Masks with Co-occurrence"
    ax.set_title(plot_title, fontsize=title_fontsize)
    
    # Extract individual masks (assuming single time step)
    mcs_mask = ds['mcs_mask'].squeeze()
    ar_mask = ds['ar_mask'].squeeze()
    etc_mask = ds['etc_mask'].squeeze()
    tc_mask = ds['tc_mask'].squeeze()
    
    # Convert individual masks to binary
    mcs_binary = xr.where(mcs_mask > 0, 1, 0)
    ar_binary = xr.where(ar_mask > 0, 1, 0)
    etc_binary = xr.where(etc_mask > 0, 1, 0)
    tc_binary = xr.where(tc_mask > 0, 1, 0)
    
    # ===== PLOT INDIVIDUAL FEATURES (BASE LAYER) =====
    if np.any(mcs_binary > 0):
        egh.healpix_show(
            mcs_binary.where(mcs_binary > 0), 
            ax=ax, 
            cmap=mpl.colors.ListedColormap([base_colors['mcs']]), 
            alpha=base_alphas['mcs'],
            vmin=0.5, vmax=1.5
        )
    
    if np.any(ar_binary > 0):
        egh.healpix_show(
            ar_binary.where(ar_binary > 0), 
            ax=ax, 
            cmap=mpl.colors.ListedColormap([base_colors['ar']]), 
            alpha=base_alphas['ar'],
            vmin=0.5, vmax=1.5
        )
    
    if np.any(etc_binary > 0):
        egh.healpix_show(
            etc_binary.where(etc_binary > 0), 
            ax=ax, 
            cmap=mpl.colors.ListedColormap([base_colors['etc']]), 
            alpha=base_alphas['etc'],
            vmin=0.5, vmax=1.5
        )
    
    if np.any(tc_binary > 0):
        egh.healpix_show(
            tc_binary.where(tc_binary > 0), 
            ax=ax, 
            cmap=mpl.colors.ListedColormap([base_colors['tc']]), 
            alpha=base_alphas['tc'],
            vmin=0.5, vmax=1.5
        )
    
    # ===== CREATE AND PLOT CO-OCCURRING FEATURE MASKS =====
    
    # 1. MCS-AR 2-way
    mcs_ar_union = xr.where(
        (ds['mcs_ar_overlap_mask'].squeeze() > 0) | (ds['ar_mcs_overlap_mask'].squeeze() > 0),
        1, 0
    )
    if np.any(mcs_ar_union > 0):
        egh.healpix_show(
            mcs_ar_union.where(mcs_ar_union > 0),
            ax=ax,
            cmap=mpl.colors.ListedColormap([cooccur_colors['mcs_ar']]),
            alpha=cooccur_alphas['mcs_ar'],
            vmin=0.5, vmax=1.5
        )
    
    # 2. MCS-ETC 2-way
    mcs_etc_union = xr.where(
        (ds['mcs_etc_overlap_mask'].squeeze() > 0) | (ds['etc_mcs_overlap_mask'].squeeze() > 0),
        1, 0
    )
    if np.any(mcs_etc_union > 0):
        egh.healpix_show(
            mcs_etc_union.where(mcs_etc_union > 0),
            ax=ax,
            cmap=mpl.colors.ListedColormap([cooccur_colors['mcs_etc']]),
            alpha=cooccur_alphas['mcs_etc'],
            vmin=0.5, vmax=1.5
        )
    
    # 3. AR-ETC 2-way
    ar_etc_union = xr.where(
        (ds['ar_etc_overlap_mask'].squeeze() > 0) | (ds['etc_ar_overlap_mask'].squeeze() > 0),
        1, 0
    )
    if np.any(ar_etc_union > 0):
        egh.healpix_show(
            ar_etc_union.where(ar_etc_union > 0),
            ax=ax,
            cmap=mpl.colors.ListedColormap([cooccur_colors['ar_etc']]),
            alpha=cooccur_alphas['ar_etc'],
            vmin=0.5, vmax=1.5
        )
    
    # 4. MCS-AR-ETC 3-way
    mcs_ar_etc_union = xr.where(
        (ds['mcs_ar_etc_overlap_mask'].squeeze() > 0) | 
        (ds['ar_mcs_etc_overlap_mask'].squeeze() > 0) | 
        (ds['etc_mcs_ar_overlap_mask'].squeeze() > 0),
        1, 0
    )
    if np.any(mcs_ar_etc_union > 0):
        egh.healpix_show(
            mcs_ar_etc_union.where(mcs_ar_etc_union > 0),
            ax=ax,
            cmap=mpl.colors.ListedColormap([cooccur_colors['mcs_ar_etc']]),
            alpha=cooccur_alphas['mcs_ar_etc'],
            vmin=0.5, vmax=1.5
        )
    
    # ===== ADD TRACK LABELS (only for co-occurring features) =====
    if show_track_labels:
        # Helper function for circular longitude averaging
        def circular_lon_mean(lon_values):
            """Calculate circular mean for longitude values to handle dateline crossing."""
            lon_rad = np.deg2rad(lon_values)
            mean_complex = np.mean(np.exp(1j * lon_rad))
            mean_lon = np.rad2deg(np.angle(mean_complex))
            if mean_lon < 0:
                mean_lon += 360
            return mean_lon
        
        # Create a combined co-occurring mask (any co-occurrence)
        any_cooccur = xr.where(
            (mcs_ar_union > 0) | (mcs_etc_union > 0) | (ar_etc_union > 0) | (mcs_ar_etc_union > 0),
            1, 0
        )
        
        # Add MCS track labels (blue) - smaller font
        mcs_tracks = np.unique(mcs_mask.values[mcs_mask.values > 0])
        mcs_fontsize = max(6, track_fontsize - 3)
        for track_id in mcs_tracks:
            track_pixels = (mcs_mask == track_id)
            cooccur_pixels = track_pixels & (any_cooccur > 0)
            if cooccur_pixels.any():
                lon_values = ds.lon.where(cooccur_pixels).values
                lat_values = ds.lat.where(cooccur_pixels).values
                valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
                if np.any(valid_mask):
                    lon_clean = lon_values[valid_mask]
                    lat_clean = lat_values[valid_mask]
                    avg_lon = circular_lon_mean(lon_clean)
                    avg_lat = np.mean(lat_clean)
                    ax.text(avg_lon, avg_lat, str(int(track_id)), 
                           transform=ccrs.PlateCarree(),
                           fontsize=mcs_fontsize, fontweight='regular', 
                           color='darkblue', ha='center', va='center')
        
        # Add ETC track labels (green)
        etc_tracks = np.unique(etc_mask.values[etc_mask.values > 0])
        for track_id in etc_tracks:
            track_pixels = (etc_mask == track_id)
            cooccur_pixels = track_pixels & (any_cooccur > 0)
            if cooccur_pixels.any():
                lon_values = ds.lon.where(cooccur_pixels).values
                lat_values = ds.lat.where(cooccur_pixels).values
                valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
                if np.any(valid_mask):
                    lon_clean = lon_values[valid_mask]
                    lat_clean = lat_values[valid_mask]
                    avg_lon = circular_lon_mean(lon_clean)
                    avg_lat = np.mean(lat_clean)
                    ax.text(avg_lon, avg_lat, str(int(track_id)), 
                           transform=ccrs.PlateCarree(),
                           fontsize=track_fontsize, fontweight='bold', 
                           color='green', ha='center', va='center')
        
        # Add AR track labels (darkorange)
        ar_tracks = np.unique(ar_mask.values[ar_mask.values > 0])
        for track_id in ar_tracks:
            track_pixels = (ar_mask == track_id)
            cooccur_pixels = track_pixels & (any_cooccur > 0)
            if cooccur_pixels.any():
                lon_values = ds.lon.where(cooccur_pixels).values
                lat_values = ds.lat.where(cooccur_pixels).values
                valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
                if np.any(valid_mask):
                    lon_clean = lon_values[valid_mask]
                    lat_clean = lat_values[valid_mask]
                    avg_lon = circular_lon_mean(lon_clean)
                    avg_lat = np.mean(lat_clean)
                    ax.text(avg_lon, avg_lat, str(int(track_id)), 
                           transform=ccrs.PlateCarree(),
                           fontsize=track_fontsize, fontweight='bold', 
                           color='darkorange', ha='center', va='center')
    
    # ===== ADD LEGENDS =====
    # Individual features legend (upper right)
    legend_base = [
        Patch(facecolor=base_colors['mcs'], alpha=base_alphas['mcs'], label='MCS'),
        Patch(facecolor=base_colors['ar'], alpha=base_alphas['ar'], label='AR'),
        Patch(facecolor=base_colors['etc'], alpha=base_alphas['etc'], label='ETC'),
        Patch(facecolor=base_colors['tc'], alpha=base_alphas['tc'], label='TC')
    ]
    
    # Co-occurring features legend (upper left)
    legend_cooccur = [
        Patch(facecolor=cooccur_colors['mcs_ar'], alpha=cooccur_alphas['mcs_ar'], label='MCS-AR'),
        Patch(facecolor=cooccur_colors['mcs_etc'], alpha=cooccur_alphas['mcs_etc'], label='MCS-ETC'),
        Patch(facecolor=cooccur_colors['ar_etc'], alpha=cooccur_alphas['ar_etc'], label='AR-ETC'),
        Patch(facecolor=cooccur_colors['mcs_ar_etc'], alpha=cooccur_alphas['mcs_ar_etc'], label='MCS-AR-ETC')
    ]
    
    # Position legends in white space from Robinson projection
    leg1 = fig.legend(handles=legend_base, 
                      bbox_to_anchor=(0.9, 0.93), 
                      loc='upper right',
                      bbox_transform=fig.transFigure,
                      fontsize=legend_fontsize, 
                      framealpha=0.9, fancybox=True, shadow=True)
    
    fig.add_artist(leg1)
    fig.legend(handles=legend_cooccur, 
               bbox_to_anchor=(0.1, 0.93), 
               loc='upper left',
               bbox_transform=fig.transFigure,
               fontsize=legend_fontsize, 
               framealpha=0.9, fancybox=True, shadow=True)
    
    # Thread-safe figure output
    if figname is not None:
        canvas = FigureCanvas(fig)
        canvas.print_png(figname)
        fig.savefig(figname, dpi=dpi, bbox_inches='tight')
        print(f"Figure saved as {figname}")
    
    plt.close(fig)
    return figname

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot feature masks for a time range with optional parallelization"
    )
    parser.add_argument("--source", required=True, 
                       help="Source name (e.g., scream)")
    parser.add_argument("--start", required=True,
                       help="Start date in YYYY-MM-DD format")
    parser.add_argument("--end", required=True,
                       help="End date in YYYY-MM-DD format")
    parser.add_argument("--parallel", type=int, default=1,
                       help="Run in parallel (0: serial, 1: parallel)")
    parser.add_argument("--figdir", default="/global/cfs/cdirs/m1867/zfeng/hk25/quicklooks/",
                       help="Output directory for figures")
    parser.add_argument("--figsize", nargs=2, type=float, default=[12, 6.75],
                       help="Figure size (width height) in inches")
    parser.add_argument("--dpi", type=int, default=200,
                       help="DPI for output figures")
    parser.add_argument("--workers", type=int, default=4,
                       help="Number of Dask workers for parallel processing")
    parser.add_argument("--plot-freq", type=str, default=None,
                       help="Override time frequency for plotting (e.g., '1H', '3H', '6H')")
    parser.add_argument("--plot-type", type=str, default="6panel", choices=["1panel", "6panel"],
                       help="Plot type: '1panel' (single panel with co-occurrence) or '6panel' (6-panel horizontal)")
    
    return parser.parse_args()

def plot_6panel_horizontal(ds, figsize=(20, 8), figname=None, dpi=200, suptitle=""):
    """
    Create a comprehensive 6-panel visualization (2x3 grid):
    - Row 0: All original features, MCS-AR overlaps, AR-ETC overlaps
    - Row 1: Isolated features, MCS-ETC overlaps, 3-way overlaps
    
    This version works with overlap masks from zarr files instead of pair lists.
    """
    
    # print("🎨 Creating comprehensive 6-panel horizontal visualization...")
    
    # Helper function for circular longitude averaging
    def circular_lon_mean(lon_values):
        """Calculate circular mean for longitude values to handle dateline crossing."""
        lon_rad = np.deg2rad(lon_values)
        mean_complex = np.mean(np.exp(1j * lon_rad))
        mean_lon = np.rad2deg(np.angle(mean_complex))
        if mean_lon < 0:
            mean_lon += 360
        return mean_lon
    
    track_fontsize = 10
    mcs_fontsize = 7
    title_fontsize = 12
    legend_fontsize = 12

    # Create GridSpec layout: 2 rows, 3 columns for equal-sized panels
    fig = plt.figure(figsize=figsize, dpi=dpi)
    gs = gridspec.GridSpec(2, 3, height_ratios=[1, 1], width_ratios=[1, 1, 1], 
                          hspace=0.14, wspace=0.05, top=0.85, bottom=0.05, left=0.05, right=0.95)
    
    # Figure suptitle
    fig.suptitle(suptitle, fontsize=title_fontsize*1.3, fontweight='bold', y=0.95)
    
    # Create projection for all axes
    projection = ccrs.Robinson(central_longitude=-135)
    
    # Define all panel positions
    ax_overview = fig.add_subplot(gs[0, 0], projection=projection)    # All original features
    ax_mcs_ar = fig.add_subplot(gs[0, 1], projection=projection)      # MCS-AR overlaps
    ax_ar_etc = fig.add_subplot(gs[0, 2], projection=projection)      # AR-ETC overlaps
    ax_isolated = fig.add_subplot(gs[1, 0], projection=projection)    # Isolated features
    ax_mcs_etc = fig.add_subplot(gs[1, 1], projection=projection)     # MCS-ETC overlaps
    ax_3way = fig.add_subplot(gs[1, 2], projection=projection)        # 3-way overlaps
    
    # Configure all axes
    all_axes = [ax_overview, ax_mcs_ar, ax_ar_etc, ax_isolated, ax_mcs_etc, ax_3way]
    for ax in all_axes:
        ax.set_global()
        ax.add_feature(cf.COASTLINE, linewidth=0.5)
    
    # ==== ROW 0, COL 0: ALL ORIGINAL FEATURES ====
    # print("  📍 Creating overview panel with all features...")
    
    # Create binary masks for visualization
    mcs_binary = (ds.mcs_mask > 0).astype(float)
    ar_binary = (ds.ar_mask > 0).astype(float)
    etc_binary = (ds.etc_mask > 0).astype(float)
    tc_binary = (ds.tc_mask > 0).astype(float)
    
    # Define colors and alphas
    overview_colors = {'mcs': 'lightblue', 'ar': 'orange', 'etc': 'limegreen', 'tc': 'red'}
    overview_alphas = {'mcs': 0.6, 'ar': 0.4, 'etc': 0.3, 'tc': 0.5}
    
    # Plot features in order (MCS first, then overlays)
    if np.any(mcs_binary > 0):
        egh.healpix_show(mcs_binary.where(mcs_binary > 0), ax=ax_overview,
                        cmap=mpl.colors.ListedColormap([overview_colors['mcs']]),
                        alpha=overview_alphas['mcs'], vmin=0.5, vmax=1.5)
    
    if np.any(ar_binary > 0):
        egh.healpix_show(ar_binary.where(ar_binary > 0), ax=ax_overview,
                        cmap=mpl.colors.ListedColormap([overview_colors['ar']]),
                        alpha=overview_alphas['ar'], vmin=0.5, vmax=1.5)
    
    if np.any(etc_binary > 0):
        egh.healpix_show(etc_binary.where(etc_binary > 0), ax=ax_overview,
                        cmap=mpl.colors.ListedColormap([overview_colors['etc']]),
                        alpha=overview_alphas['etc'], vmin=0.5, vmax=1.5)
    
    if np.any(tc_binary > 0):
        egh.healpix_show(tc_binary.where(tc_binary > 0), ax=ax_overview,
                        cmap=mpl.colors.ListedColormap([overview_colors['tc']]),
                        alpha=overview_alphas['tc'], vmin=0.5, vmax=1.5)
    
    # Add track labels for AR and ETC
    # AR track labels (orange)
    ar_tracks = np.unique(ds.ar_mask.values[ds.ar_mask.values > 0])
    for track_id in ar_tracks:
        track_pixels = (ds.ar_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_overview.text(avg_lon, avg_lat, str(int(track_id)),
                               transform=ccrs.PlateCarree(), fontsize=track_fontsize,
                               fontweight='bold', color='darkorange', ha='center', va='center')
    
    # ETC track labels (green)
    etc_tracks = np.unique(ds.etc_mask.values[ds.etc_mask.values > 0])
    for track_id in etc_tracks:
        track_pixels = (ds.etc_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_overview.text(avg_lon, avg_lat, str(int(track_id)),
                               transform=ccrs.PlateCarree(), fontsize=track_fontsize,
                               fontweight='bold', color='green', ha='center', va='center')
    
    ax_overview.set_title(f"All Original Features",
                         fontsize=title_fontsize, fontweight='bold', color='black')
    
    # ==== ROW 1, COL 0: ISOLATED FEATURES ====
    # print("  🔍 Creating isolated features panel...")
    
    # Plot isolated masks with same style as overview panel
    if np.any(ds.mcs_isolated_mask > 0):
        egh.healpix_show(ds.mcs_isolated_mask.where(ds.mcs_isolated_mask > 0), ax=ax_isolated,
                        cmap=mpl.colors.ListedColormap([overview_colors['mcs']]),
                        alpha=overview_alphas['mcs'], vmin=0.5, vmax=1.5)
    
    if np.any(ds.ar_isolated_mask > 0):
        egh.healpix_show(ds.ar_isolated_mask.where(ds.ar_isolated_mask > 0), ax=ax_isolated,
                        cmap=mpl.colors.ListedColormap([overview_colors['ar']]),
                        alpha=overview_alphas['ar'], vmin=0.5, vmax=1.5)
    
    if np.any(ds.etc_isolated_mask > 0):
        egh.healpix_show(ds.etc_isolated_mask.where(ds.etc_isolated_mask > 0), ax=ax_isolated,
                        cmap=mpl.colors.ListedColormap([overview_colors['etc']]),
                        alpha=overview_alphas['etc'], vmin=0.5, vmax=1.5)
    
    # Add track labels for isolated AR and ETC
    # AR isolated track labels (orange)
    ar_isolated_tracks = np.unique(ds.ar_isolated_mask.values[ds.ar_isolated_mask.values > 0])
    for track_id in ar_isolated_tracks:
        track_pixels = (ds.ar_isolated_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_isolated.text(avg_lon, avg_lat, str(int(track_id)),
                               transform=ccrs.PlateCarree(), fontsize=track_fontsize,
                               fontweight='bold', color='darkorange', ha='center', va='center')
    
    # ETC isolated track labels (green)
    etc_isolated_tracks = np.unique(ds.etc_isolated_mask.values[ds.etc_isolated_mask.values > 0])
    for track_id in etc_isolated_tracks:
        track_pixels = (ds.etc_isolated_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_isolated.text(avg_lon, avg_lat, str(int(track_id)),
                               transform=ccrs.PlateCarree(), fontsize=track_fontsize,
                               fontweight='bold', color='green', ha='center', va='center')
    
    # Count isolated features
    n_mcs_isolated = len(np.unique(ds.mcs_isolated_mask.values[ds.mcs_isolated_mask.values > 0]))
    n_ar_isolated = len(ar_isolated_tracks)
    n_etc_isolated = len(etc_isolated_tracks)
    
    ax_isolated.set_title(f"Isolated Features\nMCS: {n_mcs_isolated}, AR: {n_ar_isolated}, ETC: {n_etc_isolated}",
                         fontsize=title_fontsize, fontweight='bold', color='purple')
    
    # ==== ROW 0, COL 1: MCS-AR 2-WAY OVERLAPS ====
    # print("  📊 Creating MCS-AR overlap panel...")
    
    mcs_ar_overlap_tracks = np.unique(ds.mcs_ar_overlap_mask.values[ds.mcs_ar_overlap_mask.values > 0])
    ar_mcs_overlap_tracks = np.unique(ds.ar_mcs_overlap_mask.values[ds.ar_mcs_overlap_mask.values > 0])
    
    # Plot MCS-AR overlap masks
    if len(mcs_ar_overlap_tracks) > 0:
        egh.healpix_show(ds.mcs_ar_overlap_mask.where(ds.mcs_ar_overlap_mask > 0), ax=ax_mcs_ar,
                        cmap=mpl.colors.ListedColormap([overview_colors['mcs']]), alpha=0.8)
    if len(ar_mcs_overlap_tracks) > 0:
        egh.healpix_show(ds.ar_mcs_overlap_mask.where(ds.ar_mcs_overlap_mask > 0), ax=ax_mcs_ar,
                        cmap=mpl.colors.ListedColormap([overview_colors['ar']]), alpha=0.5)
    
    # Add track labels for MCS-AR overlaps
    for track_id in mcs_ar_overlap_tracks:
        track_pixels = (ds.mcs_ar_overlap_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_mcs_ar.text(avg_lon, avg_lat, str(int(track_id)), transform=ccrs.PlateCarree(),
                             fontsize=mcs_fontsize, fontweight='bold', color='darkblue', ha='center', va='center')
    
    for track_id in ar_mcs_overlap_tracks:
        track_pixels = (ds.ar_mcs_overlap_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_mcs_ar.text(avg_lon, avg_lat, str(int(track_id)), transform=ccrs.PlateCarree(),
                             fontsize=track_fontsize, fontweight='bold', color='darkorange', ha='center', va='center')
    
    ax_mcs_ar.set_title(f"MCS-AR 2-Way\nMCS: {len(mcs_ar_overlap_tracks)}, AR: {len(ar_mcs_overlap_tracks)}",
                       fontsize=title_fontsize, fontweight='bold', color='blue')
    
    # ==== ROW 0, COL 2: AR-ETC 2-WAY OVERLAPS ====
    # print("  📊 Creating AR-ETC overlap panel...")
    
    ar_etc_overlap_tracks = np.unique(ds.ar_etc_overlap_mask.values[ds.ar_etc_overlap_mask.values > 0])
    etc_ar_overlap_tracks = np.unique(ds.etc_ar_overlap_mask.values[ds.etc_ar_overlap_mask.values > 0])
    
    # Plot AR-ETC overlap masks
    if len(ar_etc_overlap_tracks) > 0:
        egh.healpix_show(ds.ar_etc_overlap_mask.where(ds.ar_etc_overlap_mask > 0), ax=ax_ar_etc,
                        cmap=mpl.colors.ListedColormap([overview_colors['ar']]), alpha=0.5)
    if len(etc_ar_overlap_tracks) > 0:
        egh.healpix_show(ds.etc_ar_overlap_mask.where(ds.etc_ar_overlap_mask > 0), ax=ax_ar_etc,
                        cmap=mpl.colors.ListedColormap([overview_colors['etc']]), alpha=0.3)

    # Add track labels for AR-ETC overlaps
    for track_id in ar_etc_overlap_tracks:
        track_pixels = (ds.ar_etc_overlap_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_ar_etc.text(avg_lon, avg_lat, str(int(track_id)), transform=ccrs.PlateCarree(),
                             fontsize=track_fontsize, fontweight='bold', color='darkorange', ha='center', va='center')
    
    for track_id in etc_ar_overlap_tracks:
        track_pixels = (ds.etc_ar_overlap_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_ar_etc.text(avg_lon, avg_lat, str(int(track_id)), transform=ccrs.PlateCarree(),
                             fontsize=track_fontsize, fontweight='bold', color='green', ha='center', va='center')
    
    ax_ar_etc.set_title(f"AR-ETC 2-Way\nAR: {len(ar_etc_overlap_tracks)}, ETC: {len(etc_ar_overlap_tracks)}",
                       fontsize=title_fontsize, fontweight='bold', color='darkorange')
    
    # ==== ROW 1, COL 1: MCS-ETC 2-WAY OVERLAPS ====
    # print("  📊 Creating MCS-ETC overlap panel...")
    
    mcs_etc_overlap_tracks = np.unique(ds.mcs_etc_overlap_mask.values[ds.mcs_etc_overlap_mask.values > 0])
    etc_mcs_overlap_tracks = np.unique(ds.etc_mcs_overlap_mask.values[ds.etc_mcs_overlap_mask.values > 0])
    
    # Plot MCS-ETC overlap masks
    if len(mcs_etc_overlap_tracks) > 0:
        egh.healpix_show(ds.mcs_etc_overlap_mask.where(ds.mcs_etc_overlap_mask > 0), ax=ax_mcs_etc,
                        cmap=mpl.colors.ListedColormap([overview_colors['mcs']]), alpha=0.8)
    if len(etc_mcs_overlap_tracks) > 0:
        egh.healpix_show(ds.etc_mcs_overlap_mask.where(ds.etc_mcs_overlap_mask > 0), ax=ax_mcs_etc,
                        cmap=mpl.colors.ListedColormap([overview_colors['etc']]), alpha=0.3)

    # Add track labels for MCS-ETC overlaps
    for track_id in mcs_etc_overlap_tracks:
        track_pixels = (ds.mcs_etc_overlap_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_mcs_etc.text(avg_lon, avg_lat, str(int(track_id)), transform=ccrs.PlateCarree(),
                               fontsize=mcs_fontsize, fontweight='bold', color='darkblue', ha='center', va='center')
    
    for track_id in etc_mcs_overlap_tracks:
        track_pixels = (ds.etc_mcs_overlap_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_mcs_etc.text(avg_lon, avg_lat, str(int(track_id)), transform=ccrs.PlateCarree(),
                               fontsize=track_fontsize, fontweight='bold', color='green', ha='center', va='center')
    
    ax_mcs_etc.set_title(f"MCS-ETC 2-Way\nMCS: {len(mcs_etc_overlap_tracks)}, ETC: {len(etc_mcs_overlap_tracks)}",
                        fontsize=title_fontsize, fontweight='bold', color='limegreen')
    
    # ==== ROW 1, COL 2: 3-WAY OVERLAPS ====
    # print("  📊 Creating 3-way overlap panel...")
    
    mcs_3way_overlap_tracks = np.unique(ds.mcs_ar_etc_overlap_mask.values[ds.mcs_ar_etc_overlap_mask.values > 0])
    ar_3way_overlap_tracks = np.unique(ds.ar_mcs_etc_overlap_mask.values[ds.ar_mcs_etc_overlap_mask.values > 0])
    etc_3way_overlap_tracks = np.unique(ds.etc_mcs_ar_overlap_mask.values[ds.etc_mcs_ar_overlap_mask.values > 0])
    
    # Plot 3-way overlap masks
    if len(mcs_3way_overlap_tracks) > 0:
        egh.healpix_show(ds.mcs_ar_etc_overlap_mask.where(ds.mcs_ar_etc_overlap_mask > 0), ax=ax_3way,
                        cmap=mpl.colors.ListedColormap([overview_colors['mcs']]), alpha=0.8)
    if len(ar_3way_overlap_tracks) > 0:
        egh.healpix_show(ds.ar_mcs_etc_overlap_mask.where(ds.ar_mcs_etc_overlap_mask > 0), ax=ax_3way,
                        cmap=mpl.colors.ListedColormap([overview_colors['ar']]), alpha=0.5)
    if len(etc_3way_overlap_tracks) > 0:
        egh.healpix_show(ds.etc_mcs_ar_overlap_mask.where(ds.etc_mcs_ar_overlap_mask > 0), ax=ax_3way,
                        cmap=mpl.colors.ListedColormap([overview_colors['etc']]), alpha=0.3)
    
    # Add track labels for 3-way overlaps
    for track_id in mcs_3way_overlap_tracks:
        track_pixels = (ds.mcs_ar_etc_overlap_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_3way.text(avg_lon, avg_lat, str(int(track_id)), transform=ccrs.PlateCarree(),
                           fontsize=mcs_fontsize, fontweight='bold', color='darkblue', ha='center', va='center')
    
    for track_id in ar_3way_overlap_tracks:
        track_pixels = (ds.ar_mcs_etc_overlap_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_3way.text(avg_lon, avg_lat, str(int(track_id)), transform=ccrs.PlateCarree(),
                           fontsize=track_fontsize, fontweight='bold', color='darkorange', ha='center', va='center')
    
    for track_id in etc_3way_overlap_tracks:
        track_pixels = (ds.etc_mcs_ar_overlap_mask == track_id)
        if track_pixels.any():
            lon_values = ds.lon.where(track_pixels).values
            lat_values = ds.lat.where(track_pixels).values
            valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
            if np.any(valid_mask):
                avg_lon = circular_lon_mean(lon_values[valid_mask])
                avg_lat = np.mean(lat_values[valid_mask])
                ax_3way.text(avg_lon, avg_lat, str(int(track_id)), transform=ccrs.PlateCarree(),
                           fontsize=track_fontsize, fontweight='bold', color='green', ha='center', va='center')
    
    ax_3way.set_title(f"3-Way Overlaps\nMCS: {len(mcs_3way_overlap_tracks)}, AR: {len(ar_3way_overlap_tracks)}, ETC: {len(etc_3way_overlap_tracks)}",
                     fontsize=title_fontsize, fontweight='bold', color='red')
    
    # Add legend
    legend_elements = [
        Patch(facecolor=overview_colors['mcs'], alpha=overview_alphas['mcs'], label='MCS'),
        Patch(facecolor=overview_colors['ar'], alpha=overview_alphas['ar'], label='AR'),
        Patch(facecolor=overview_colors['etc'], alpha=overview_alphas['etc'], label='ETC'),
        Patch(facecolor=overview_colors['tc'], alpha=overview_alphas['tc'], label='TC')
    ]
    fig.legend(handles=legend_elements, bbox_to_anchor=(0.05, 0.95), loc='upper left', bbox_transform=fig.transFigure,
               fontsize=legend_fontsize, framealpha=0.9, fancybox=False, shadow=False)

    # Thread-safe figure output
    if figname is not None:
        canvas = FigureCanvas(fig)
        canvas.print_png(figname)
        fig.savefig(figname, dpi=dpi, bbox_inches='tight')
        print(f"Figure saved as {figname}")
    
    plt.close(fig)
    return figname


def process_single_time(data_path, time_step_str, source_name, figdir, figsize, dpi, plot_type="6panel"):
    """
    Process a single time step and create the plot.
    This function is designed to work with Dask delayed execution.
    Instead of passing the entire dataset, we pass the path and load only the needed time step.
    
    Parameters:
    - data_path: string path to the zarr dataset
    - time_step_str: string representation of the time step (YYYY-MM-DDTHH)
    - source_name: source name for title and filename
    - figdir: output directory for figures
    - figsize: figure size tuple
    - dpi: output DPI
    - plot_type: '1panel' or '6panel' (default: '6panel')
    
    Returns:
    - success: int (1 for success, 0 for failure)
    """
    try:
        # Load dataset and select only the specific time step
        ds = xr.open_zarr(data_path, consolidated=True)
        
        # Select single time step
        _ds = ds.sel(time=time_step_str, method='nearest')
        _time = ds.time.sel(time=time_step_str, method='nearest')
        # Attach Healpix coordinates
        _ds = _ds.pipe(egh.attach_coords)

        # Close the full dataset to free memory
        ds.close()
        
        # Format time for title and filename
        time_label = _time.dt.strftime('%Y-%m-%d %H:%M UTC').item()
        time_str = _time.dt.strftime('%Y%m%d_%H%M').item()
        title_with_time = f"Co-occurring Features {time_label} ({source_name.upper()})"
        
        # Choose plotting function and filename based on plot_type
        if plot_type == "1panel":
            figname = f'{figdir}{source_name}_cofmasks_1panel_{time_str}.png'
            _figname = plot_features_with_cooccurrence(
                _ds, 
                title=title_with_time, 
                figsize=figsize,
                figname=figname,
                dpi=dpi,
                show_track_labels=True,
                track_fontsize=8,
            )
        else:  # 6panel
            figname = f'{figdir}{source_name}_cofmasks_6panel_{time_str}.png'
            _figname = plot_6panel_horizontal(
                _ds, 
                suptitle=title_with_time, 
                figsize=figsize,
                figname=figname,
                dpi=dpi,
            )
        
        return 1
        
    except Exception as e:
        print(f"Error processing time {time_step_str}: {e}")
        # Ensure cleanup even on error
        plt.close('all')
        return 0

def main():
    """Main function to orchestrate the plotting process."""
    # Configure matplotlib for memory efficiency
    plt.ioff()  # Turn off interactive mode
    
    # Parse command line arguments
    args = parse_args()
    
    # Set up paths and parameters
    source_name = args.source
    start_date = args.start
    end_date = args.end
    run_parallel = args.parallel
    figdir = args.figdir
    figsize = tuple(args.figsize)
    dpi = args.dpi
    n_workers = args.workers
    plot_type = args.plot_type
    
    # Create output directory
    os.makedirs(figdir, exist_ok=True)
    
    # Set up data paths
    root_dir = "/pscratch/sd/w/wcmca1/hackathon/cof_masks/"
    in_dir = f"{root_dir}/{source_name}_cofmasks_hp8_v1.zarr"
    
    print(f"Loading data from: {in_dir}")
    print(f"Plot type: {plot_type}")
    
    # Load dataset only to get time information and validate
    try:
        ds = xr.open_zarr(in_dir, consolidated=True)
        # Get available times to validate date range
        available_times = ds.time.values
        print(f"  ✅ Dataset loaded successfully")
        print(f"  Time steps: {len(ds.time)}")
        print(f"  Data variables: {list(ds.data_vars)}")
        print(f"  Spatial dimensions: {dict(ds.dims)}")
        ds.close()  # Close immediately to free memory
    except Exception as e:
        print(f"  ❌ Error loading dataset: {e}")
        return
    
    # Determine time frequency for plotting
    if args.plot_freq is not None:
        # Use user-provided frequency
        freq_str = args.plot_freq
        print(f"  Using user-specified frequency: {freq_str}")
    elif len(available_times) > 1:
        # Calculate average time frequency from the dataset
        time_diffs = pd.to_datetime(available_times[1:]) - pd.to_datetime(available_times[:-1])
        avg_freq = time_diffs.mean()
        
        # Convert to frequency string (e.g., '3H', '1H', '6H')
        total_seconds = avg_freq.total_seconds()
        if total_seconds >= 3600:  # >= 1 hour
            hours = int(total_seconds / 3600)
            freq_str = f'{hours}H'
        elif total_seconds >= 60:  # >= 1 minute
            minutes = int(total_seconds / 60)
            freq_str = f'{minutes}T'  # T is pandas notation for minutes
        else:
            seconds = int(total_seconds)
            freq_str = f'{seconds}S'
        
        print(f"  Calculated time frequency from dataset: {freq_str} (avg interval: {avg_freq})")
    else:
        freq_str = '6H'  # fallback to default
        print(f"  Using default frequency: {freq_str} (insufficient data points to calculate)")
    
    # Create date range using determined frequency
    time_range = pd.date_range(start=start_date, end=end_date, freq=freq_str)
    
    print(f"Creating plots from {start_date} to {end_date}")
    print(f"Total time steps: {len(time_range)}")
    print(f"Figure size: {figsize}, DPI: {dpi}")
    print(f"Output directory: {figdir}")
    
    if run_parallel == 0:
        # Serial processing
        print("Running in serial mode...")
        success_count = 0
        for i, time_step in enumerate(time_range):
            tm_str = time_step.strftime('%Y-%m-%dT%H')
            result = process_single_time(in_dir, tm_str, source_name, figdir, figsize, dpi, plot_type)
            success_count += result
            
            # Progress update
            if (i + 1) % 10 == 0 or i == len(time_range) - 1:
                print(f"Completed {i + 1}/{len(time_range)} plots ({success_count} successful)")
                
    else:
        # Parallel processing with Dask
        print(f"Running in parallel mode with {n_workers} workers...")
        
        # Set up Dask cluster
        client = setup_dask_client(parallel=True, n_workers=n_workers, threads_per_worker=1)
        print(f"Dask cluster ready with {n_workers} workers")
        
        try:
            # Process tasks in batches to avoid large graph serialization
            print("Starting parallel computation...")
            batch_size = min(n_workers * 2, 32)  # Process in batches of 2x workers or 32, whichever is smaller
            total_time_steps = len(time_range)
            success_count = 0
            completed_count = 0
            
            from dask.distributed import as_completed
            
            # Process time steps in batches
            for batch_start in range(0, total_time_steps, batch_size):
                batch_end = min(batch_start + batch_size, total_time_steps)
                batch_times = time_range[batch_start:batch_end]
                
                print(f"Processing batch {batch_start//batch_size + 1}/{-(-total_time_steps//batch_size)}: "
                      f"time steps {batch_start+1}-{batch_end}")
                
                # Submit batch of tasks
                futures = []
                for time_step in batch_times:
                    tm_str = time_step.strftime('%Y-%m-%dT%H')
                    future = client.submit(
                        process_single_time, 
                        in_dir, tm_str, source_name, figdir, figsize, dpi, plot_type
                    )
                    futures.append(future)
                
                # Process batch results
                for future in as_completed(futures):
                    try:
                        result = future.result()
                        success_count += result
                        completed_count += 1
                        
                        # Progress update every 10 plots
                        if completed_count % 10 == 0:
                            print(f"Progress: {completed_count}/{total_time_steps} plots completed ({success_count} successful)")
                            
                    except Exception as e:
                        print(f"Task failed: {e}")
                        completed_count += 1
            
            print(f"Completed all {total_time_steps} time steps ({success_count} successful)")
            
        finally:
            # Clean up Dask resources
            if client is not None:
                client.close()
    
    print("Finished creating all plots!")

if __name__ == "__main__":
    main()
