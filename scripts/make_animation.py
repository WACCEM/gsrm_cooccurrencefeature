#!/usr/bin/env python
"""
Make multi-feature tracking mask animation

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import os
import subprocess
import pandas as pd
import tempfile

###############################################################################################
# Script parameters
###############################################################################################

# Choose feature type: "all" or "cof"
feature_type = "cof"  # Change this to "all" or "cof"

# For "cof" type, choose plot style: "1panel" or "6panel"
cof_plot_type = "1panel"  # Only used when feature_type="cof"

# Script parameters
# source_name = "scream"
source_name = "IMERGv7"
# source_name = "icon_d3hp003"

# start_date = "2020-02-25T00"
# end_date = "2020-03-01T23"
# start_date = "2020-12-01T00"
# end_date = "2020-12-31T23"
start_date = "2019-01-01T00"
end_date = "2019-01-31T23"

parallel_mode = 1
n_workers = 64
output_dpi = 200

# Optional: Override time frequency for plotting (set to None to auto-calculate)
# Examples: '1h', '3h', '6h', '12h'
plot_freq = '6h'  # Set to desired frequency string or None for auto-calculation

# Execution control options
run_plotting = True   # Set to False to skip plotting and use existing PNG files
run_ffmpeg = True     # Set to False to skip animation creation (plotting only)

# FFmpeg animation parameters
input_framerate = 1    # (frames per second) - how fast to transition between frames
output_framerate = 10  # (frames per second) - video playback speed (lower values = smaller file size, e.g., 24 fps is cinema)
video_quality = 20     # CRF value (lower = better quality, range 0-51, 18-28 is good)
output_width = None    # Optional: output video width in pixels (e.g., 1920, 1280). None = use original size. Height scales automatically to maintain aspect ratio.

# Set paths and options based on feature type
if feature_type == "all":
    # All masks configuration
    fig_basename = "allmasks"
    figdir = f"/global/cfs/cdirs/m1867/zfeng/hk25/quicklooks/{source_name}/"
    plotting_code = "plot_feature_masks.py"
    fig_width = 14
    fig_height = 8
elif feature_type == "cof":
    # Co-occurrence masks configuration
    if cof_plot_type == "1panel":
        fig_basename = "cofmasks_1panel"
        fig_width = 12
        fig_height = 6
    else:  # 6panel
        fig_basename = "cofmasks_6panel"
        fig_width = 20
        fig_height = 8
    figdir = f"/global/cfs/cdirs/m1867/zfeng/hk25/quicklooks_cof/{source_name}/"
    plotting_code = "plot_cooccurrence_masks.py"

# Animation parameters
animation_dir = "/global/cfs/cdirs/m1867/zfeng/hk25/animations/"
start_date_str = start_date.split('T')[0]  # Extract YYYY-MM-DD
end_date_str = end_date.split('T')[0]      # Extract YYYY-MM-DD
animation_filename = f"{animation_dir}{source_name}_{fig_basename}_{start_date_str}_{end_date_str}.mp4"

###############################################################################################
# Main execution
###############################################################################################

print("Make multi-feature tracking mask animation")
print(f"Feature type: {feature_type}")
print(f"Source: {source_name}")
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
    cmd = [
        'python', plotting_code,
        '--source', source_name,
        '--start', start_date,
        '--end', end_date,
        '--parallel', str(parallel_mode),
        '--workers', str(n_workers),
        '--figdir', figdir,
        '--figsize', str(fig_width), str(fig_height),
        '--dpi', str(output_dpi)
    ]

    # Add plot frequency if specified
    if plot_freq is not None:
        cmd.extend(['--plot-freq', plot_freq])
        print(f"Using custom plot frequency: {plot_freq}")
    else:
        print("Using auto-calculated plot frequency from dataset")
    
    # Add plot type for co-occurrence plots
    if feature_type == "cof":
        cmd.extend(['--plot-type', cof_plot_type])
        print(f"Co-occurrence plot type: {cof_plot_type}")

    print(f"Command: {' '.join(cmd)}")
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
        freq_str = '6H'  # Default fallback used in plotting scripts

    time_range = pd.date_range(start=start_date, end=end_date, freq=freq_str)

    # Create list of expected PNG filenames
    expected_files = []
    for time_step in time_range:
        time_str = time_step.strftime('%Y%m%d_%H%M')
        png_filename = f'{figdir}{source_name}_{fig_basename}_{time_str}.png'
        
        # Check if file actually exists
        if os.path.exists(png_filename):
            expected_files.append(png_filename)

    print(f"Found {len(expected_files)} PNG files within date range")
    print(f"Time range: {start_date} to {end_date} (freq: {freq_str})")

    if len(expected_files) == 0:
        print("❌ No PNG files found within the specified date range!")
        print(f"   Check directory: {figdir}")
        print(f"   Expected pattern: {source_name}_{fig_basename}_YYYYMMDD_HHMM.png")
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