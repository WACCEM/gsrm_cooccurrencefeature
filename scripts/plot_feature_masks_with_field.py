#!/usr/bin/env python
"""
Plot feature masks with background field (e.g., IVT, precipitation) and optional contours 
for a time range with Dask parallelization.

Usage:
python plot_feature_masks_with_field.py --source scream --catalog-model scream_ne120_inst \
    --start 2019-10-31 --end 2019-11-02 [options]

Optional arguments:
--parallel 0 (serial), 1 (parallel) - default: 1
--catalog-url URL - default: https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml
--current-location LOCATION - default: NERSC
--catalog-params PARAMS - default: {"zoom": 8}
--figdir output_directory - default: /global/cfs/cdirs/m1867/zfeng/hk25/quicklooks/
--figsize width height - figure size in inches, default: 12 6
--dpi DPI for output figures, default: 200
--workers number of dask workers, default: 4
--plot-freq time frequency for plotting (e.g., '1H', '3H', '6H')
--plot-mcs plot MCS masks (default: False)
--no-ar disable AR masks (enabled by default)
--no-etc disable ETC masks (enabled by default)
--no-tc disable TC masks (enabled by default)
--show-mcs-labels show MCS track labels (default: False)
--no-ar-labels disable AR track labels (enabled by default)
--no-etc-labels disable ETC track labels (enabled by default)
--no-tc-labels disable TC track labels (enabled by default)
--shade-var variable name for shading (default: IVT)
--contour-var variable name for contours (default: psl)
--no-shading disable background shading
--no-contour disable contours

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import argparse
import os
import json
import cftime
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib as mpl
# Set non-GUI backend for thread safety
mpl.use('agg')
import matplotlib.pyplot as plt
# Configure matplotlib to be more memory-efficient
plt.rcParams['figure.max_open_warning'] = 0
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
import intake
import colormaps as cmaps

# Suppress warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def get_colormap(cmap_name):
    """Get colormap from name, supporting matplotlib and cmaps."""
    if cmap_name is None:
        return None
    
    # Try cmaps module first
    if hasattr(cmaps, cmap_name):
        return getattr(cmaps, cmap_name)
    
    # Handle reversed colormaps
    if cmap_name.endswith('_r'):
        base_name = cmap_name[:-2]
        if hasattr(cmaps, base_name):
            return getattr(cmaps, base_name).reversed()
        # Try matplotlib
        try:
            return plt.cm.get_cmap(cmap_name)
        except:
            pass
    
    # Try matplotlib colormaps
    try:
        return plt.cm.get_cmap(cmap_name)
    except:
        print(f"Warning: Colormap '{cmap_name}' not found, using default")
        return None

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot feature masks with background field for a time range with optional parallelization"
    )
    parser.add_argument("--source", required=True, 
                       help="Source name (e.g., scream)")
    parser.add_argument("--catalog-model", required=True,
                       help="Catalog model name for field data (e.g., scream_ne120_inst)")
    parser.add_argument("--start", required=True,
                       help="Start date in YYYY-MM-DD format")
    parser.add_argument("--end", required=True,
                       help="End date in YYYY-MM-DD format")
    parser.add_argument("--parallel", type=int, default=1,
                       help="Run in parallel (0: serial, 1: parallel)")
    parser.add_argument("--catalog-url", 
                       default="https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml",
                       help="URL to the intake catalog")
    parser.add_argument("--current-location", default="NERSC",
                       help="Current location for catalog access")
    parser.add_argument("--catalog-params", type=str, default='{"zoom": 8}',
                       help='Catalog parameters as JSON string (e.g., \'{"zoom": 8}\')')
    parser.add_argument("--figdir", default="/global/cfs/cdirs/m1867/zfeng/hk25/quicklooks/",
                       help="Output directory for figures")
    parser.add_argument("--figsize", nargs=2, type=float, default=[12, 6],
                       help="Figure size (width height) in inches")
    parser.add_argument("--dpi", type=int, default=200,
                       help="DPI for output figures")
    parser.add_argument("--workers", type=int, default=8,
                       help="Number of Dask workers for parallel processing")
    parser.add_argument("--plot-freq", type=str, default=None,
                       help="Override time frequency for plotting (e.g., '1H', '3H', '6H')")
    
    # Mask plotting flags
    parser.add_argument("--plot-mcs", action="store_true", default=False,
                       help="Plot MCS masks")
    parser.add_argument("--no-ar", action="store_false", dest="plot_ar", default=True,
                       help="Disable AR masks (enabled by default)")
    parser.add_argument("--no-etc", action="store_false", dest="plot_etc", default=True,
                       help="Disable ETC masks (enabled by default)")
    parser.add_argument("--no-tc", action="store_false", dest="plot_tc", default=True,
                       help="Disable TC masks (enabled by default)")
    
    # Track label flags
    parser.add_argument("--show-mcs-labels", action="store_true", default=False,
                       help="Show MCS track labels")
    parser.add_argument("--no-ar-labels", action="store_false", dest="show_ar_labels", default=True,
                       help="Disable AR track labels (enabled by default)")
    parser.add_argument("--no-etc-labels", action="store_false", dest="show_etc_labels", default=True,
                       help="Disable ETC track labels (enabled by default)")
    parser.add_argument("--no-tc-labels", action="store_false", dest="show_tc_labels", default=True,
                       help="Disable TC track labels (enabled by default)")
    
    # Field plotting options
    parser.add_argument("--shade-var", type=str, default="IVT",
                       help="Variable name for background shading (must exist in field dataset)")
    parser.add_argument("--shade-cmap", type=str, default=None,
                       help="Colormap name for shading (e.g., 'viridis', 'gray_r', 'blue_medb717b')")
    parser.add_argument("--shade-levels", type=str, default=None,
                       help="Shade levels as 'min,max,step' (e.g., '80,320,5' for np.arange(80, 320, 5))")
    parser.add_argument("--contour-var", type=str, default="psl",
                       help="Variable name for contours (must exist in field dataset)")
    parser.add_argument("--no-shading", action="store_true", default=False,
                       help="Disable background shading")
    parser.add_argument("--no-contour", action="store_true", default=False,
                       help="Disable contours")
    parser.add_argument("--mcs-cmap", type=str, default=None,
                       help="Colormap name for MCS masks (e.g., 'tab20', 'cet_g_bw_minc_minl')")
    
    return parser.parse_args()


def plot_feature_mask_with_field(
    ds, dsm, title="", figsize=(12, 6), 
    base_colors=None, base_alphas=None, text_colors=None,
    figname=None, dpi=200, 
    title_fontsize=14, legend_fontsize=12,
    # Mask plotting flags
    plot_mcs=False, plot_ar=True, plot_etc=True, plot_tc=True,
    # Track label flags
    show_mcs_labels=False, show_ar_labels=True, show_etc_labels=True, show_tc_labels=True,
    track_fontsize=12,
    # Shading options
    plot_shading=True, shade_var='IVT', shade_cmap=None, shade_levels=None,
    # Contour options
    plot_contour=True, contour_var='psl', contour_cmap=None, contour_levels=None, contour_lw=0.5,
    # MCS colormap and levels
    mcs_cmap=None, mcs_levels=None,
):
    """
    Plot feature masks with optional shading and contours.
    
    Parameters:
    -----------
    ds : xr.Dataset
        Dataset containing fields to plot (e.g., IVT, psl)
    dsm : xr.Dataset
        Dataset containing mask variables
    title : str
        Plot title
    figsize : tuple
        Figure size (width, height)
    base_colors : dict
        Colors for each feature type {'mcs': color, 'ar': color, 'etc': color, 'tc': color}
    base_alphas : dict
        Alpha values for each feature type and shading
    text_colors : dict
        Colors for track labels
    figname : str
        Output filename (if None, figure not saved)
    dpi : int
        Figure resolution
    title_fontsize : int
        Title font size
    legend_fontsize : int
        Legend font size
    plot_mcs : bool
        Whether to plot MCS masks
    plot_ar : bool
        Whether to plot AR masks
    plot_etc : bool
        Whether to plot ETC masks
    plot_tc : bool
        Whether to plot TC masks
    show_mcs_labels : bool
        Whether to show MCS track labels
    show_ar_labels : bool
        Whether to show AR track labels
    show_etc_labels : bool
        Whether to show ETC track labels
    show_tc_labels : bool
        Whether to show TC track labels
    track_fontsize : int
        Font size for track labels
    plot_shading : bool
        Whether to plot shaded background field
    shade_var : str
        Variable name for shading (must exist in ds)
    shade_cmap : colormap or None
        Colormap for shading
    shade_levels : array-like or None
        Levels for shading normalization
    plot_contour : bool
        Whether to plot contours
    contour_var : str
        Variable name for contours (must exist in ds)
    contour_cmap : str or None
        Colormap for contours
    contour_levels : array-like or None
        Contour levels
    contour_lw : float
        Contour line width
    mcs_cmap : colormap or None
        Colormap specifically for MCS masks
    mcs_levels : array-like or None
        Track number levels for MCS masks to ensure consistent coloring across time steps
        
    Returns:
    --------
    figname : str or None
        Output filename
    """
    
    # Set default colors
    if base_colors is None:
        base_colors = {
            'mcs': 'lightskyblue',
            'ar': 'darkorange', 
            'etc': 'limegreen',
            'tc': 'red',
        }

    if text_colors is None:
        text_colors = {
            'mcs': 'navy',
            'ar': 'darkorange', 
            'etc': 'green',
            'tc': 'maroon',
        }
    
    if base_alphas is None:
        base_alphas = {
            'mcs': 0.8,
            'ar': 0.5,
            'etc': 0.5,
            'tc': 0.6,
            'shade': 0.9,
        }
    
    # Set default MCS colormap
    if mcs_cmap is None:
        mcs_cmap = cmaps.cet_g_bw_minc_minl
    
    # Set default shading colormap and levels
    if shade_cmap is None:
        shade_cmap = cmaps.blue_medb717b
    if shade_levels is None:
        shade_levels = np.arange(0, 1201, 10)
    
    # Set default contour levels and colormap
    if contour_levels is None:
        contour_levels = np.arange(95000, 105001, 1000)
    if contour_cmap is None:
        contour_cmap = 'plasma'

    # Extract individual masks (assuming single time step)
    mcs_mask = dsm['mcs_mask'].squeeze()
    ar_mask = dsm['ar_mask'].squeeze()
    etc_mask = dsm['etc_mask'].squeeze()
    tc_mask = dsm['tc_mask'].squeeze()

    # Convert individual masks to binary
    ar_binary = xr.where(ar_mask > 0, 1, 0)
    etc_binary = xr.where(etc_mask > 0, 1, 0)
    tc_binary = xr.where(tc_mask > 0, 1, 0)

    # Create base map
    projection = ccrs.Robinson(central_longitude=-135)
    fig = plt.figure(figsize=figsize, dpi=dpi, facecolor='w')
    ax = plt.subplot(111, projection=projection)
    ax.set_global()
    ax.add_feature(cf.COASTLINE, linewidth=0.8)

    # Gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 10}
    gl.ylabel_style = {"size": 10}

    # Set title
    plot_title = title if title else ""
    ax.set_title(plot_title, fontsize=title_fontsize)

    # Plot shaded background field if requested
    if plot_shading:
        shade_field = ds[shade_var].squeeze()
        norm = mpl.colors.BoundaryNorm(boundaries=shade_levels, ncolors=shade_cmap.N)
        alpha_shade = base_alphas.get('shade', 0.9)
        img_shade = egh.healpix_show(shade_field, cmap=shade_cmap, norm=norm, ax=ax, alpha=alpha_shade, zorder=1)

    # Plot masks based on flags
    zorder = 2
    
    # Plot ETC mask
    if plot_etc and np.any(etc_binary > 0):
        img_etc = egh.healpix_show(etc_binary.where(etc_binary > 0), ax=ax, 
                    cmap=mpl.colors.ListedColormap([base_colors['etc']]), 
                    alpha=base_alphas['etc'],
                    vmin=0.5, vmax=1.5, zorder=zorder)
        zorder += 1

    # Plot AR mask
    if plot_ar and np.any(ar_binary > 0):
        img_ar = egh.healpix_show(ar_binary.where(ar_binary > 0), ax=ax, 
                    cmap=mpl.colors.ListedColormap([base_colors['ar']]), 
                    alpha=base_alphas['ar'],
                    vmin=0.5, vmax=1.5, zorder=zorder)
        zorder += 1

    # Plot MCS mask
    if plot_mcs and np.any(mcs_mask > 0):
        # Use BoundaryNorm if mcs_levels provided for consistent track coloring
        if mcs_levels is not None and len(mcs_levels) > 0:
            # Handle case where there are more tracks than colormap colors
            if len(mcs_levels) > 255:
                # print(f"Warning: {len(mcs_levels)} MCS tracks exceed typical colormap size (255)")
                # Use modulo to cycle through colors for large track counts
                mcs_mask_cycled = xr.where(mcs_mask > 0, (mcs_mask - 1) % 255 + 1, 0)
                norm = mpl.colors.BoundaryNorm(boundaries=np.arange(0, 257), ncolors=mcs_cmap.N)
                img_mcs = egh.healpix_show(mcs_mask_cycled.where(mcs_mask_cycled > 0), ax=ax, 
                            cmap=mcs_cmap, norm=norm,
                            alpha=base_alphas['mcs'], zorder=zorder)
            else:
                norm = mpl.colors.BoundaryNorm(boundaries=mcs_levels, ncolors=mcs_cmap.N)
                img_mcs = egh.healpix_show(mcs_mask.where(mcs_mask > 0), ax=ax, 
                            cmap=mcs_cmap, norm=norm,
                            alpha=base_alphas['mcs'], zorder=zorder)
        else:
            img_mcs = egh.healpix_show(mcs_mask.where(mcs_mask > 0), ax=ax, 
                        cmap=mcs_cmap, 
                        alpha=base_alphas['mcs'], zorder=zorder)
        zorder += 1

    # Plot TC mask
    if plot_tc and np.any(tc_binary > 0):
        img_tc = egh.healpix_show(tc_binary.where(tc_binary > 0), ax=ax, 
                    cmap=mpl.colors.ListedColormap([base_colors['tc']]), 
                    alpha=base_alphas['tc'],
                    vmin=0.5, vmax=1.5, zorder=zorder)
        zorder += 1

    # Plot contours if requested
    if plot_contour:
        contour_field = ds[contour_var].squeeze()
        img_cont = egh.healpix_contour(contour_field, cmap=contour_cmap, levels=contour_levels, 
                                       linewidths=contour_lw, ax=ax, zorder=zorder)
        zorder += 1

    # Helper function for circular longitude averaging (handles dateline crossing)
    def circular_lon_mean(lon_values):
        """Calculate circular mean for longitude values to handle dateline crossing."""
        lon_rad = np.deg2rad(lon_values)
        mean_complex = np.mean(np.exp(1j * lon_rad))
        mean_lon = np.rad2deg(np.angle(mean_complex))
        if mean_lon < 0:
            mean_lon += 360
        return mean_lon
    
    # Add track labels for features if requested
    label_zorder = zorder
    
    # Add ETC track labels
    if show_etc_labels:
        etc_tracks = np.unique(etc_mask.values[etc_mask.values > 0])
        for track_id in etc_tracks:
            track_pixels = (etc_mask == track_id)
            if track_pixels.any():
                lon_values = ds.lon.where(track_pixels).values
                lat_values = ds.lat.where(track_pixels).values
                
                valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
                if np.any(valid_mask):
                    lon_clean = lon_values[valid_mask]
                    lat_clean = lat_values[valid_mask]
                    
                    avg_lon = circular_lon_mean(lon_clean)
                    avg_lat = np.mean(lat_clean)
                    
                    ax.text(avg_lon, avg_lat, str(int(track_id)), 
                           transform=ccrs.PlateCarree(), zorder=label_zorder,
                           fontsize=track_fontsize, fontweight='bold', 
                           color=text_colors['etc'], ha='center', va='center')
    
    # Add AR track labels
    if show_ar_labels:
        ar_tracks = np.unique(ar_mask.values[ar_mask.values > 0])
        for track_id in ar_tracks:
            track_pixels = (ar_mask == track_id)
            if track_pixels.any():
                lon_values = ds.lon.where(track_pixels).values
                lat_values = ds.lat.where(track_pixels).values
                
                valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
                if np.any(valid_mask):
                    lon_clean = lon_values[valid_mask]
                    lat_clean = lat_values[valid_mask]
                    
                    avg_lon = circular_lon_mean(lon_clean)
                    avg_lat = np.mean(lat_clean)
                    
                    ax.text(avg_lon, avg_lat, str(int(track_id)), 
                           transform=ccrs.PlateCarree(), zorder=label_zorder,
                           fontsize=track_fontsize, fontweight='bold', 
                           color=text_colors['ar'], ha='center', va='center')
    
    # Add TC track labels
    if show_tc_labels:
        tc_tracks = np.unique(tc_mask.values[tc_mask.values > 0])
        for track_id in tc_tracks:
            track_pixels = (tc_mask == track_id)
            if track_pixels.any():
                lon_values = ds.lon.where(track_pixels).values
                lat_values = ds.lat.where(track_pixels).values
                
                valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
                if np.any(valid_mask):
                    lon_clean = lon_values[valid_mask]
                    lat_clean = lat_values[valid_mask]
                    
                    avg_lon = circular_lon_mean(lon_clean)
                    avg_lat = np.mean(lat_clean)
                    
                    ax.text(avg_lon, avg_lat, str(int(track_id)), 
                           transform=ccrs.PlateCarree(), zorder=label_zorder,
                           fontsize=track_fontsize, fontweight='bold', 
                           color=text_colors['tc'], ha='center', va='center')
    
    # Add MCS track labels
    if show_mcs_labels:
        mcs_tracks = np.unique(mcs_mask.values[mcs_mask.values > 0])
        for track_id in mcs_tracks:
            track_pixels = (mcs_mask == track_id)
            if track_pixels.any():
                lon_values = ds.lon.where(track_pixels).values
                lat_values = ds.lat.where(track_pixels).values
                
                valid_mask = ~np.isnan(lon_values) & ~np.isnan(lat_values)
                if np.any(valid_mask):
                    lon_clean = lon_values[valid_mask]
                    lat_clean = lat_values[valid_mask]
                    
                    avg_lon = circular_lon_mean(lon_clean)
                    avg_lat = np.mean(lat_clean)
                    
                    ax.text(avg_lon, avg_lat, str(int(track_id)), 
                           transform=ccrs.PlateCarree(), zorder=label_zorder,
                           fontsize=track_fontsize, fontweight='bold', 
                           color=text_colors['mcs'], ha='center', va='center')
    
    # Build legend dynamically based on what's being plotted
    legend_elements = []
    if plot_mcs:
        legend_elements.append(Patch(facecolor=base_colors['mcs'], alpha=base_alphas['mcs'], label='MCS'))
    if plot_ar:
        legend_elements.append(Patch(facecolor=base_colors['ar'], alpha=base_alphas['ar'], label='AR'))
    if plot_etc:
        legend_elements.append(Patch(facecolor=base_colors['etc'], alpha=base_alphas['etc'], label='ETC'))
    if plot_tc:
        legend_elements.append(Patch(facecolor=base_colors['tc'], alpha=base_alphas['tc'], label='TC'))
    
    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper right', fontsize=legend_fontsize, 
                  framealpha=0.9, fancybox=True, shadow=True)
    
    # Thread-safe figure output
    if figname is not None:
        canvas = FigureCanvas(fig)
        canvas.print_png(figname)
        fig.savefig(figname, dpi=dpi, bbox_inches='tight')
        print(f"Figure saved as {figname}")
    
    plt.close(fig)
    return figname


def process_single_time(mask_path, field_info, time_step_str, source_name, 
                       figdir, figsize, dpi, shade_var, contour_var,
                       plot_mcs, plot_ar, plot_etc, plot_tc,
                       show_mcs_labels, show_ar_labels, show_etc_labels, show_tc_labels,
                       plot_shading, plot_contour,
                       shade_cmap=None, shade_levels=None, mcs_cmap=None, mcs_levels=None):
    """
    Process a single time step and create the plot.
    This function is designed to work with Dask delayed execution.
    
    Parameters:
    - mask_path: string path to the zarr mask dataset
    - field_info: dict with catalog info {'catalog_url', 'current_location', 'catalog_model', 'catalog_params'}
    - time_step_str: string representation of the time step (YYYY-MM-DDTHH)
    - source_name: source name for title and filename
    - figdir: output directory for figures
    - figsize: figure size tuple
    - dpi: output DPI
    - shade_var: variable name for shading
    - contour_var: variable name for contours
    - plot_mcs, plot_ar, plot_etc, plot_tc: mask plotting flags
    - show_mcs_labels, show_ar_labels, show_etc_labels, show_tc_labels: label flags
    - plot_shading: whether to plot shading
    - plot_contour: whether to plot contours
    
    Returns:
    - success: int (1 for success, 0 for failure)
    """
    try:
        # Parse time string to datetime object
        time_obj = pd.to_datetime(time_step_str)
        
        # Load mask dataset and select only the specific time step
        dsm = xr.open_zarr(mask_path, consolidated=True)
        
        # Handle cftime calendars by converting time_obj if needed
        if hasattr(dsm.time.values[0], 'calendar'):
            # Dataset uses cftime, convert datetime to cftime
            calendar_type = dsm.time.values[0].calendar
            time_to_select = cftime.datetime(
                time_obj.year, time_obj.month, time_obj.day, 
                time_obj.hour, time_obj.minute, time_obj.second,
                calendar=calendar_type
            )
        else:
            time_to_select = time_obj
        
        _dsm = dsm.sel(time=time_to_select, method='nearest')
        _time = dsm.time.sel(time=time_to_select, method='nearest')
        dsm.close()
        
        # Load field data from catalog for this specific time step
        cat = intake.open_catalog(field_info['catalog_url'])[field_info['current_location']]
        field_dataset = cat[field_info['catalog_model']](**field_info['catalog_params']).to_dask().pipe(egh.attach_coords, signed_lon=True)
        
        # Handle cftime calendars for field dataset
        if hasattr(field_dataset.time.values[0], 'calendar'):
            calendar_type = field_dataset.time.values[0].calendar
            field_time_to_select = cftime.datetime(
                time_obj.year, time_obj.month, time_obj.day, 
                time_obj.hour, time_obj.minute, time_obj.second,
                calendar=calendar_type
            )
        else:
            field_time_to_select = time_obj
        
        _ds = field_dataset.sel(time=field_time_to_select, method='nearest')
        
        # Compute IVT if shade_var is 'IVT' and it doesn't exist
        if shade_var == 'IVT' and 'IVT' not in _ds:
            if 'uivt' in _ds and 'vivt' in _ds:
                _uivt = _ds.uivt
                _vivt = _ds.vivt
                _ivt = np.sqrt(_uivt**2 + _vivt**2)
                _ds['IVT'] = _ivt
            else:
                print(f"Warning: IVT requested but uivt/vivt not found for {time_step_str}")
                return 0
        
        # Attach Healpix coordinates to both datasets
        _ds = _ds.pipe(egh.attach_coords)
        _dsm = _dsm.pipe(egh.attach_coords)
        
        # Format time for title and filename
        time_label = _time.dt.strftime('%Y-%m-%d %H:%M UTC').item()
        time_str = _time.dt.strftime('%Y%m%d_%H%M').item()
        title_with_time = f"{source_name.upper()} {time_label}"
        figname = f'{figdir}{source_name}_masks_{shade_var}_{time_str}.png'
        
        # Create plot
        _figname = plot_feature_mask_with_field(
            _ds, _dsm,
            title=title_with_time,
            figsize=figsize,
            figname=figname,
            dpi=dpi,
            title_fontsize=18,
            legend_fontsize=12,
            track_fontsize=10,
            plot_mcs=plot_mcs,
            plot_ar=plot_ar,
            plot_etc=plot_tc,
            plot_tc=plot_tc,
            show_mcs_labels=show_mcs_labels,
            show_ar_labels=show_ar_labels,
            show_etc_labels=show_etc_labels,
            show_tc_labels=show_tc_labels,
            plot_shading=plot_shading,
            shade_var=shade_var,
            shade_cmap=shade_cmap,
            shade_levels=shade_levels,
            plot_contour=plot_contour,
            contour_var=contour_var,
            mcs_cmap=mcs_cmap,
            mcs_levels=mcs_levels,
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
    catalog_model = args.catalog_model
    start_date = args.start
    end_date = args.end
    run_parallel = args.parallel
    figdir = args.figdir
    figsize = tuple(args.figsize)
    dpi = args.dpi
    n_workers = args.workers
    
    # Parse catalog parameters
    try:
        catalog_params = json.loads(args.catalog_params)
    except json.JSONDecodeError as e:
        print(f"Error parsing catalog-params JSON: {e}")
        return
    
    # Parse colormaps
    shade_cmap = get_colormap(args.shade_cmap) if args.shade_cmap else None
    mcs_cmap = get_colormap(args.mcs_cmap) if args.mcs_cmap else None
    
    # Parse shade levels
    shade_levels = None
    if args.shade_levels:
        try:
            parts = args.shade_levels.split(',')
            if len(parts) == 3:
                vmin, vmax, step = map(float, parts)
                shade_levels = np.arange(vmin, vmax, step)
            else:
                print(f"Warning: shade-levels must be 'min,max,step', got '{args.shade_levels}'")
        except ValueError as e:
            print(f"Warning: Error parsing shade-levels: {e}")
    
    # Create output directory
    os.makedirs(figdir, exist_ok=True)
    
    # Set up data paths
    root_dir = "/pscratch/sd/w/wcmca1/hackathon/cof_masks/"
    mask_path = f"{root_dir}/{source_name}_cofmasks_hp8_v1.zarr"
    
    print(f"Loading mask data from: {mask_path}")
    
    # Load mask dataset only to get time information and validate
    try:
        dsm = xr.open_zarr(mask_path, consolidated=True)
        available_times = dsm.time.values
        print(f"  ✅ Mask dataset loaded successfully")
        print(f"  Time steps: {len(dsm.time)}")
        print(f"  Data variables: {list(dsm.data_vars)}")

        # If MCS plotting is requested, get the min/max MCS track numbers within the time range
        mcs_levels = None
        if args.plot_mcs:
            # Select only the requested time range for computing min/max
            dsm_subset = dsm.sel(time=slice(start_date, end_date))
            mcs_min = int(dsm_subset['mcs_mask'].min().values.item())
            mcs_max = int(dsm_subset['mcs_mask'].max().values.item())
            mcs_levels = np.arange(mcs_min, mcs_max + 2)  # +2 for boundary norm (min to max+1)
            print(f"  MCS tracks range from {mcs_min} to {mcs_max} ({len(mcs_levels)-1} unique tracks) within requested period")
            dsm_subset.close()
        dsm.close()
    except Exception as e:
        print(f"  ❌ Error loading mask dataset: {e}")
        return
    
    # Prepare field data catalog information (to be passed to workers)
    print(f"Field data catalog configuration:")
    print(f"  Catalog URL: {args.catalog_url}")
    print(f"  Location: {args.current_location}")
    print(f"  Model: {catalog_model}")
    print(f"  Parameters: {catalog_params}")
    
    # Create field_info dict for passing to workers (lightweight, no dask dataset)
    field_info = {
        'catalog_url': args.catalog_url,
        'current_location': args.current_location,
        'catalog_model': catalog_model,
        'catalog_params': catalog_params
    }
    
    # Validate catalog access by loading one time step
    try:
        print(f"  Validating catalog access...")
        cat = intake.open_catalog(args.catalog_url)[args.current_location]
        ds_field_test = cat[catalog_model](**catalog_params).to_dask().pipe(egh.attach_coords, signed_lon=True)
        print(f"  ✅ Field dataset catalog validated successfully")
        print(f"  Time steps: {len(ds_field_test.time)}")
        print(f"  Data variables: {list(ds_field_test.data_vars)}")
        del ds_field_test  # Free memory
    except Exception as e:
        print(f"  ❌ Error accessing field dataset catalog: {e}")
        return
    
    # Determine time frequency for plotting
    if args.plot_freq is not None:
        freq_str = args.plot_freq
        print(f"  Using user-specified frequency: {freq_str}")
    elif len(available_times) > 1:
        time_diffs = pd.to_datetime(available_times[1:]) - pd.to_datetime(available_times[:-1])
        avg_freq = time_diffs.mean()
        
        total_seconds = avg_freq.total_seconds()
        if total_seconds >= 3600:
            hours = int(total_seconds / 3600)
            freq_str = f'{hours}h'
        elif total_seconds >= 60:
            minutes = int(total_seconds / 60)
            freq_str = f'{minutes}min'
        else:
            seconds = int(total_seconds)
            freq_str = f'{seconds}s'
        
        print(f"  Calculated time frequency from dataset: {freq_str} (avg interval: {avg_freq})")
    else:
        freq_str = '3h'
        print(f"  Using default frequency: {freq_str}")
    
    # Create date range using determined frequency
    time_range = pd.date_range(start=start_date, end=end_date, freq=freq_str)
    
    print(f"Creating plots from {start_date} to {end_date}")
    print(f"Total time steps: {len(time_range)}")
    print(f"Figure size: {figsize}, DPI: {dpi}")
    print(f"Output directory: {figdir}")
    print(f"Shading: {not args.no_shading} ({args.shade_var if not args.no_shading else 'N/A'})")
    print(f"Contours: {not args.no_contour} ({args.contour_var if not args.no_contour else 'N/A'})")
    
    if run_parallel == 0:
        # Serial processing
        print("Running in serial mode...")
        success_count = 0
        for i, time_step in enumerate(time_range):
            tm_str = time_step.strftime('%Y-%m-%dT%H')
            result = process_single_time(
                mask_path, field_info, tm_str, source_name, 
                figdir, figsize, dpi,
                args.shade_var, args.contour_var,
                args.plot_mcs, args.plot_ar, args.plot_etc, args.plot_tc,
                args.show_mcs_labels, args.show_ar_labels, args.show_etc_labels, args.show_tc_labels,
                not args.no_shading, not args.no_contour,
                shade_cmap, shade_levels, mcs_cmap, mcs_levels
            )
            success_count += result
            
            if (i + 1) % 10 == 0 or i == len(time_range) - 1:
                print(f"Completed {i + 1}/{len(time_range)} plots ({success_count} successful)")
                
    else:
        # Parallel processing with Dask
        print(f"Running in parallel mode with {n_workers} workers...")
        
        client = setup_dask_client(parallel=True, n_workers=n_workers, threads_per_worker=1)
        print(f"Dask cluster ready with {n_workers} workers")
        
        try:
            print("Starting parallel computation...")
            batch_size = min(n_workers * 2, 32)
            total_time_steps = len(time_range)
            success_count = 0
            completed_count = 0
            
            from dask.distributed import as_completed
            
            for batch_start in range(0, total_time_steps, batch_size):
                batch_end = min(batch_start + batch_size, total_time_steps)
                batch_times = time_range[batch_start:batch_end]
                
                print(f"Processing batch {batch_start//batch_size + 1}/{-(-total_time_steps//batch_size)}: "
                      f"time steps {batch_start+1}-{batch_end}")
                
                futures = []
                for time_step in batch_times:
                    tm_str = time_step.strftime('%Y-%m-%dT%H')
                    future = client.submit(
                        process_single_time,
                        mask_path, field_info, tm_str, source_name,
                        figdir, figsize, dpi,
                        args.shade_var, args.contour_var,
                        args.plot_mcs, args.plot_ar, args.plot_etc, args.plot_tc,
                        args.show_mcs_labels, args.show_ar_labels, args.show_etc_labels, args.show_tc_labels,
                        not args.no_shading, not args.no_contour,
                        shade_cmap, shade_levels, mcs_cmap, mcs_levels
                    )
                    futures.append(future)
                
                for future in as_completed(futures):
                    try:
                        result = future.result()
                        success_count += result
                        completed_count += 1
                        
                        if completed_count % 10 == 0:
                            print(f"Progress: {completed_count}/{total_time_steps} plots completed ({success_count} successful)")
                            
                    except Exception as e:
                        print(f"Task failed: {e}")
                        completed_count += 1
            
            print(f"Completed all {total_time_steps} time steps ({success_count} successful)")
            
        finally:
            if client is not None:
                client.close()
    
    print("Finished creating all plots!")


if __name__ == "__main__":
    main()
