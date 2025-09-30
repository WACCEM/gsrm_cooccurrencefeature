#!/bin/bash
###############################################################################################
# Make multi-feature tracking mask animation
###############################################################################################

# Script parameters
# source_name="scream"
source_name="IMERGv7"
start_date="2019-11-24T00"
end_date="2019-11-30T23"
parallel_mode=1
n_workers=32
fig_width=14
fig_height=8
output_dpi=200

# Figure directory
# fig_basename="allmasks"
fig_basename="cofmasks"
# figdir="/global/cfs/cdirs/m1867/zfeng/hk25/quicklooks/"${source_name}"/"
figdir="/global/cfs/cdirs/m1867/zfeng/hk25/quicklooks_cof/"${source_name}"/"
# Animation parameters
animation_dir="/global/cfs/cdirs/m1867/zfeng/hk25/animations/"
# animation_filename="${animation_dir}${source_name}_allmasks_$(echo ${start_date} | cut -d'T' -f1).mp4"
animation_filename="${animation_dir}${source_name}_${fig_basename}_$(echo ${start_date} | cut -d'T' -f1).mp4"

# Activate Python environment
source activate /global/common/software/m1867/python/hackathon
# Create directories if they don't exist
mkdir -p ${animation_dir}
mkdir -p ${figdir}

# plotting_code="plot_feature_masks.py"
plotting_code="plot_cooccurrence_masks.py"

# Run plotting script
echo "Running parallel processing with ${n_workers} workers"
python ${plotting_code} \
    --source ${source_name} \
    --start ${start_date} \
    --end ${end_date} \
    --parallel ${parallel_mode} \
    --workers ${n_workers} \
    --figdir ${figdir} \
    --figsize ${fig_width} ${fig_height} \
    --dpi ${output_dpi}

# Make animation using ffmpeg
echo "Making animation from quicklook plots..."
ffmpeg -framerate 2 -pattern_type glob -i "${figdir}${source_name}_${fig_basename}_*.png" \
    -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
    -c:v libx264 -r 10 -crf 20 -pix_fmt yuv420p \
    -y ${animation_filename}
echo "View animation here: ${animation_filename}"
