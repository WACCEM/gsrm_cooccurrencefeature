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
    
    return parser.parse_args()

def plot_all_feature_masks(ds, title="", figsize=(12, 6), colors=None, figname=None, dpi=200, 
                          alphas=None, title_fontsize=14, legend_fontsize=12,
                          show_track_labels=True, track_fontsize=8):
    """
    Plot all 4 feature masks at a single time step with thread-safe matplotlib.
    
    Parameters:
    - ds: xarray Dataset containing the mask variables (assumes single time step)
    - title: optional title for the plot
    - figsize: figure size tuple
    - colors: dict with color names for each feature type
    - figname: if provided, save figure as PNG with this filename
    - dpi: DPI for figure display and saving (default 200)
    - alphas: dict with alpha values for each feature type
    - title_fontsize: font size for the plot title (default 14)
    - legend_fontsize: font size for the legend (default 12)
    - show_track_labels: whether to show track ID labels for ETC and AR features (default True)
    - track_fontsize: font size for track labels (default 8)
    
    Returns:
    - figname: filename of the saved figure
    """
    # Set default colors
    if colors is None:
        colors = {
            'mcs': 'lightskyblue',
            'ar': 'darkorange', 
            'etc': 'limegreen',
            'tc': 'crimson'
        }
    
    # Set default alphas
    if alphas is None:
        alphas = {
            'mcs': 0.8,
            'ar': 0.5,
            'etc': 0.3,
            'tc': 0.6
        }
    
    # Create base map (integrated from worldmap_base)
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
    plot_title = title if title else "Feature Masks"
    ax.set_title(plot_title, fontsize=title_fontsize)
    
    # Extract masks (assuming single time step)
    mcs_mask = ds['mcs_mask'].squeeze()
    ar_mask = ds['ar_mask'].squeeze()
    etc_mask = ds['etc_mask'].squeeze()
    tc_mask = ds['tc_mask'].squeeze()
    
    # Convert masks to binary (all features = 1, no features = 0)
    mcs_binary = xr.where(mcs_mask > 0, 1, 0)
    ar_binary = xr.where(ar_mask > 0, 1, 0)
    etc_binary = xr.where(etc_mask > 0, 1, 0)
    tc_binary = xr.where(tc_mask > 0, 1, 0)
    
    # Plot MCS mask with customizable color and alpha
    if np.any(mcs_binary > 0):
        img_mcs = egh.healpix_show(
            mcs_binary.where(mcs_binary > 0), 
            ax=ax, 
            cmap=mpl.colors.ListedColormap([colors['mcs']]), 
            alpha=alphas['mcs'],
            vmin=0.5, 
            vmax=1.5
        )
    
    # Plot AR mask with customizable color and alpha
    if np.any(ar_binary > 0):
        img_ar = egh.healpix_show(
            ar_binary.where(ar_binary > 0), 
            ax=ax, 
            cmap=mpl.colors.ListedColormap([colors['ar']]), 
            alpha=alphas['ar'],
            vmin=0.5, 
            vmax=1.5
        )
    
    # Plot ETC mask with customizable color and alpha
    if np.any(etc_binary > 0):
        img_etc = egh.healpix_show(
            etc_binary.where(etc_binary > 0), 
            ax=ax, 
            cmap=mpl.colors.ListedColormap([colors['etc']]), 
            alpha=alphas['etc'],
            vmin=0.5, 
            vmax=1.5
        )
    
    # Plot TC mask with customizable color and alpha
    if np.any(tc_binary > 0):
        img_tc = egh.healpix_show(
            tc_binary.where(tc_binary > 0), 
            ax=ax, 
            cmap=mpl.colors.ListedColormap([colors['tc']]), 
            alpha=alphas['tc'],
            vmin=0.5, 
            vmax=1.5
        )

    # Add track labels for ETC and AR features if requested
    if show_track_labels:
        # Helper function for circular longitude averaging (handles dateline crossing)
        def circular_lon_mean(lon_values):
            """Calculate circular mean for longitude values to handle dateline crossing."""
            # Convert to radians
            lon_rad = np.deg2rad(lon_values)
            # Calculate circular mean using complex numbers
            mean_complex = np.mean(np.exp(1j * lon_rad))
            # Convert back to degrees and ensure 0-360 range
            mean_lon = np.rad2deg(np.angle(mean_complex))
            if mean_lon < 0:
                mean_lon += 360
            return mean_lon
        
        # Add ETC track labels (green)
        etc_tracks = np.unique(etc_mask.values[etc_mask.values > 0])
        for track_id in etc_tracks:
            track_pixels = (etc_mask == track_id)
            if track_pixels.any():
                # Calculate average lat/lon for this track using circular statistics for longitude
                lon_values = ds.lon.where(track_pixels).values
                lat_values = ds.lat.where(track_pixels).values
                
                # Remove NaN values
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
            if track_pixels.any():
                # Calculate average lat/lon for this track using circular statistics for longitude
                lon_values = ds.lon.where(track_pixels).values
                lat_values = ds.lat.where(track_pixels).values
                
                # Remove NaN values
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
    
    # Add legend
    legend_elements = [
        Patch(facecolor=colors['mcs'], alpha=alphas['mcs'], label='MCS'),
        Patch(facecolor=colors['ar'], alpha=alphas['ar'], label='AR'),
        Patch(facecolor=colors['etc'], alpha=alphas['etc'], label='ETC'),
        Patch(facecolor=colors['tc'], alpha=alphas['tc'], label='TC')
    ]
    # fig.legend(handles=legend_elements, bbox_to_anchor=(0.99, 0.99), loc='upper right', bbox_transform=fig.transFigure,
    #            fontsize=legend_fontsize, framealpha=0.9, fancybox=True, shadow=True)
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


def process_single_time(data_path, time_step_str, source_name, figdir, figsize, dpi):
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
        title_with_time = f"{source_name.upper()} {time_label}"
        figname = f'{figdir}{source_name}_allmasks_{time_str}.png'
        
        # Set up consistent styling
        custom_colors = {
            'mcs': 'lightskyblue',
            'ar': 'darkorange', 
            'etc': 'limegreen',
            'tc': 'crimson'
        }

        custom_alphas = {
            'mcs': 0.8,    
            'ar': 0.5,     
            'etc': 0.3,    
            'tc': 0.6      
        }
        
        # Create plot
        _figname = plot_all_feature_masks(
            _ds, 
            title=title_with_time, 
            figsize=figsize,
            colors=custom_colors,
            alphas=custom_alphas,
            title_fontsize=18,
            legend_fontsize=12,
            figname=figname,
            dpi=dpi,
            show_track_labels=True,
            track_fontsize=8,
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
    
    # Create output directory
    os.makedirs(figdir, exist_ok=True)
    
    # Set up data paths
    root_dir = "/pscratch/sd/w/wcmca1/hackathon/allmasks/"
    in_dir = f"{root_dir}/{source_name}_allmasks_hp8_v1.zarr"
    
    print(f"Loading data from: {in_dir}")
    
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
        freq_str = '3H'  # fallback to default
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
            result = process_single_time(in_dir, tm_str, source_name, figdir, figsize, dpi)
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
                        in_dir, tm_str, source_name, figdir, figsize, dpi
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
