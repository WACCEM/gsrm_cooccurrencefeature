#!/usr/bin/env python
"""
Make multi-feature tracking mask animation with background field (IVT, RLUT, etc.)

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import os
import subprocess
import pandas as pd
import tempfile

###############################################################################################
# Script parameters
###############################################################################################

# Choose feature type: "etc_ar_tc" or "mcs"
# feature_type = "etc_ar_tc"
feature_type = "mcs"

# Script parameters
source_name = "scream"

# start_date = "2019-08-02T00"
# end_date = "2019-08-31T23"
start_date = "2019-08-03T00"
end_date = "2019-08-08T23"

parallel_mode = 1
n_workers = 32
output_dpi = 200

# Optional: Override time frequency for plotting (set to None to auto-calculate)
# Examples: '1h', '3h', '6h', '12h'
plot_freq = '6h'  # Set to desired frequency string or None for auto-calculation

# Execution control options
run_plotting = False    # Set to False to skip plotting and use existing PNG files
run_ffmpeg = True     # Set to False to skip animation creation (plotting only)

# FFmpeg animation parameters
input_framerate = 1    # (frames per second) - how fast to transition between frames
output_framerate = 10  # (frames per second) - video playback speed (lower values = smaller file size, e.g., 24 fps is cinema)
video_quality = 20     # CRF value (lower = better quality, range 0-51, 18-28 is good)
output_width = 1280    # Optional: output video width in pixels (e.g., 1920, 1280). None = use original size. Height scales automatically to maintain aspect ratio.

# Set paths and options based on feature type
if feature_type == "etc_ar_tc":
    # ETC + AR + TC configuration with IVT shading and PSL contours
    fig_basename = "masks_ivt_psl"
    figdir = f"/global/cfs/cdirs/m1867/zfeng/hk25/quicklooks_field/{source_name}/"
    catalog_model = "scream_ne120_inst"
    fig_width = 12
    fig_height = 6
    
    # Feature mask options
    plot_mcs = False
    plot_ar = True
    plot_etc = True
    plot_tc = True
    show_mcs_labels = False
    show_ar_labels = True
    show_etc_labels = True
    show_tc_labels = True
    
    # Field options
    shade_var = "IVT"
    shade_cmap = None  # Use default (colormaps.blue_medb717b)
    shade_levels = None  # Use default (0-1200)
    contour_var = "psl"
    no_shading = False
    no_contour = False
    mcs_cmap = None
    
elif feature_type == "mcs":
    # MCS configuration with RLUT shading
    fig_basename = "masks_mcs_rlut"
    figdir = f"/global/cfs/cdirs/m1867/zfeng/hk25/quicklooks_field/{source_name}/"
    catalog_model = "scream2D_hrly"
    fig_width = 12
    fig_height = 6
    
    # Feature mask options
    plot_mcs = True
    plot_ar = False
    plot_etc = False
    plot_tc = False
    show_mcs_labels = False
    show_ar_labels = False
    show_etc_labels = False
    show_tc_labels = True
    
    # Field options
    shade_var = "rlut"
    shade_cmap = "gray_r"
    shade_levels = "80,320,5"
    contour_var = "psl"
    no_shading = False
    no_contour = True
    mcs_cmap = "cet_g_bw_minc_minl"
else:
    print(f"❌ Error: Unknown feature_type '{feature_type}'")
    print("   Valid options: 'etc_ar_tc' or 'mcs'")
    exit(1)

# Animation parameters
animation_dir = "/global/cfs/cdirs/m1867/zfeng/hk25/animations/"
start_date_str = start_date.split('T')[0]  # Extract YYYY-MM-DD
end_date_str = end_date.split('T')[0]      # Extract YYYY-MM-DD
animation_filename = f"{animation_dir}{source_name}_{fig_basename}_{start_date_str}_{end_date_str}.mp4"

plotting_code = "plot_feature_masks_with_field.py"

###############################################################################################
# Main execution
###############################################################################################

print("Make multi-feature tracking mask animation with background field")
print(f"Feature type: {feature_type}")
print(f"Source: {source_name}")
print(f"Catalog model: {catalog_model}")
print(f"Date range: {start_date} to {end_date}")
print(f"Execution mode: Plotting={'✅' if run_plotting else '❌'}, FFmpeg={'✅' if run_ffmpeg else '❌'}")

# Create directories if they don't exist
if run_ffmpeg:
    os.makedirs(animation_dir, exist_ok=True)
if run_plotting:
    os.makedirs(figdir, exist_ok=True)

# Run plotting script
if run_plotting:
    print(f"📊 Running plotting script with {n_workers} workers...")
    
    # Build base command
    cmd = [
        'python', plotting_code,
        '--source', source_name,
        '--catalog-model', catalog_model,
        '--start', start_date,
        '--end', end_date,
        '--parallel', str(parallel_mode),
        '--workers', str(n_workers),
        '--figdir', figdir,
        '--figsize', str(fig_width), str(fig_height),
        '--dpi', str(output_dpi),
        '--shade-var', shade_var,
        '--contour-var', contour_var,
    ]

    # Add plot frequency if specified
    if plot_freq is not None:
        cmd.extend(['--plot-freq', plot_freq])
        print(f"Using custom plot frequency: {plot_freq}")
    else:
        print("Using auto-calculated plot frequency from dataset")
    
    # Add mask plotting flags
    if plot_mcs:
        cmd.append('--plot-mcs')
    if not plot_ar:
        cmd.append('--no-ar')
    if not plot_etc:
        cmd.append('--no-etc')
    if not plot_tc:
        cmd.append('--no-tc')
    
    # Add track label flags
    if show_mcs_labels:
        cmd.append('--show-mcs-labels')
    if not show_ar_labels:
        cmd.append('--no-ar-labels')
    if not show_etc_labels:
        cmd.append('--no-etc-labels')
    if not show_tc_labels:
        cmd.append('--no-tc-labels')
    
    # Add field options
    if no_shading:
        cmd.append('--no-shading')
    if no_contour:
        cmd.append('--no-contour')
    
    if shade_cmap is not None:
        cmd.extend(['--shade-cmap', shade_cmap])
    if shade_levels is not None:
        cmd.extend(['--shade-levels', shade_levels])
    if mcs_cmap is not None:
        cmd.extend(['--mcs-cmap', mcs_cmap])

    print(f"Command: {' '.join(cmd)}")
    print(f"Configuration:")
    print(f"  Masks: MCS={plot_mcs}, AR={plot_ar}, ETC={plot_etc}, TC={plot_tc}")
    print(f"  Labels: MCS={show_mcs_labels}, AR={show_ar_labels}, ETC={show_etc_labels}, TC={show_tc_labels}")
    print(f"  Shading: {shade_var} {'(disabled)' if no_shading else ''}")
    if shade_cmap:
        print(f"    Colormap: {shade_cmap}")
    if shade_levels:
        print(f"    Levels: {shade_levels}")
    print(f"  Contours: {contour_var} {'(disabled)' if no_contour else ''}")
    if mcs_cmap:
        print(f"  MCS colormap: {mcs_cmap}")
    
    result = subprocess.run(cmd)

    if result.returncode != 0:
        print(f"❌ Error: Plotting script failed with exit code {result.returncode}")
        exit(1)
        
    print("✅ Plotting completed successfully!")
else:
    print("⏭️  Skipping plotting - using existing PNG files")

# Make animation using ffmpeg
if run_ffmpeg:
    print("🎬 Creating animation from PNG files...")

    # Generate expected time range based on the same logic as plotting
    if plot_freq is not None:
        freq_str = plot_freq
    else:
        freq_str = '6h'  # Default fallback used in plotting scripts

    time_range = pd.date_range(start=start_date, end=end_date, freq=freq_str)

    # Create list of expected PNG filenames
    # Note: The plotting script uses pattern: {source}_masks_{shade_var}_{time_str}.png
    expected_files = []
    for time_step in time_range:
        time_str = time_step.strftime('%Y%m%d_%H%M')
        png_filename = f'{figdir}{source_name}_masks_{shade_var}_{time_str}.png'
        
        # Check if file actually exists
        if os.path.exists(png_filename):
            expected_files.append(png_filename)

    print(f"Found {len(expected_files)} PNG files within date range")
    print(f"Time range: {start_date} to {end_date} (freq: {freq_str})")

    if len(expected_files) == 0:
        print("❌ No PNG files found within the specified date range!")
        print(f"   Check directory: {figdir}")
        print(f"   Expected pattern: {source_name}_masks_{shade_var}_YYYYMMDD_HHMM.png")
        exit(1)

    # Create a temporary file list for FFmpeg
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        temp_filelist = f.name
        for png_file in sorted(expected_files):
            f.write(f"file '{png_file}'\n")

    try:
        # Build scale filter based on output_width parameter
        if output_width is not None:
            # Scale to specified width, maintain aspect ratio with height divisible by 2
            scale_filter = f'scale={output_width}:-2'
        else:
            # Just ensure dimensions are even (required for H.264 encoding)
            scale_filter = 'scale=trunc(iw/2)*2:trunc(ih/2)*2'
        
        # Use FFmpeg concat demuxer with file list
        ffmpeg_cmd = [
            'ffmpeg',
            '-f', 'concat',
            '-safe', '0',
            '-r', str(input_framerate),  # Input framerate
            '-i', temp_filelist,
            '-vf', scale_filter,
            '-c:v', 'libx264',
            '-r', str(output_framerate),  # Output framerate
            '-crf', str(video_quality),
            '-pix_fmt', 'yuv420p',
            '-y', animation_filename
        ]
        
        print(f"🎬 Animation settings:")
        print(f"   Input framerate: {input_framerate} fps (PNG reading speed)")
        print(f"   Output framerate: {output_framerate} fps (video playback speed)")
        print(f"   Video quality (CRF): {video_quality} (lower=better)")
        print(f"   Output width: {output_width if output_width else 'original'} pixels")
        print(f"FFmpeg command: {' '.join(ffmpeg_cmd)}")
        print(f"Using {len(expected_files)} PNG files from {expected_files[0]} to {expected_files[-1]}")
        
        result = subprocess.run(ffmpeg_cmd)
        
    finally:
        # Clean up temporary file
        os.unlink(temp_filelist)

    if result.returncode == 0:
        print(f"✅ Animation created successfully!")
        print(f"🎬 View animation here: {animation_filename}")
    else:
        print(f"❌ Error: FFmpeg failed with exit code {result.returncode}")
        
else:
    print("⏭️  Skipping animation creation - PNG files ready for manual processing")
