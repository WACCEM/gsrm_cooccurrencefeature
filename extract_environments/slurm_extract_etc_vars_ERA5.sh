#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
##SBATCH -q shared
##SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH -J era5
#SBATCH -A m1867
#SBATCH --array=0-2%10
#SBATCH --output=logs/extract_etc_era5_%A_%a.log
#SBATCH --mail-user=zhe.feng@pnnl.gov
#SBATCH --mail-type=FAIL,END

# ===== JOB ARRAY STRUCTURE FOR ERA5 =====
# This script uses SLURM job arrays for ERA5 ETC extraction
# ERA5 tracks are on structured lat/lon grid (not unstructured HEALPix mesh)
# --array=0-2%10 means:
#   - 3 tasks total (one per variable, indices 0-2)
#   - %10 limits to 10 simultaneous tasks
# Benefits:
#   - Parallel processing: Process multiple variables simultaneously
#   - Fault tolerance: If one variable fails, others continue
#   - Resource isolation: Each variable gets dedicated resources
#   - Easy restart: Can resubmit only failed array indices
# 
# The concurrent setting can be overridden at job submission:
#   sbatch --array=0-2%5 slurm_extract_etc_vars_ERA5.sh
# To restart failed tasks: 
#   sbatch --array=1 slurm_extract_etc_vars_ERA5.sh
# ================================

# module load python
# module list
source activate /global/common/software/m1867/python/hackathon

# Create logs directory if it doesn't exist
mkdir -p logs

# Set up paths and parameters
ROOT_DIR="/pscratch/sd/w/wcmca1/hackathon/etc_tracks"
TRACK_FILE="${ROOT_DIR}/era5.etc_stitched_nodes.txt"
OUTPUT_DIR="/pscratch/sd/w/wcmca1/hackathon/etc_data/era5/single_vars/"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# ===== PARAMETERS TO CUSTOMIZE =====
# Model and catalog settings
CATALOG_URL="/global/homes/f/feng045/program/hackathon/catalog/NERSC/main.yaml"
CURRENT_LOCATION=""  # Empty for local catalog
CATALOG_MODEL="era5_3h"  # ERA5 model (3-hourly data)
CATALOG_PARAMS='{"zoom": 8}'

# ===== VARIABLE CONFIGURATION =====
# 
# 2D VARIABLES (no pressure dimension):
#   pr ps psl tas d2m u10 v10 sstk prw
#   mcs_ar_etc_overlap_mask, etc_mcs_ar_overlap_mask, ar_mcs_etc_overlap_mask
#
# 3D VARIABLES (require --pressure_levels):
#   ua, va, omega, hus (specific humidity)
#   Note: For 3D variables, specify pressure levels below
#
# NOTE: For ERA5 + 'pr' variable, IMERG V7 precipitation will be used automatically
#
# Set variables to extract (this corresponds to SLURM array indices)
# E.g., for 3 variables, do: --array=0-2
VARIABLES=(
    # "ua" "va" "omega" "hur" "zg"
#   "pr" "ps" "psl" "tas" "d2m" "u10" "v10" "sstk" "prw"
#   "ua" "va" "omega" "hus" "zg"  # 3D variables (need PRESSURE_LEVELS)
  "mcs_ar_etc_overlap_mask"  # COF mask (2D, requires --cof_mask flag)
  "ar_mcs_etc_overlap_mask"  # COF mask (2D, requires --cof_mask flag)
  "etc_mcs_ar_overlap_mask"  # COF mask (2D, requires --cof_mask flag)
)

# 3D variable options (for pressure level data)
# Leave empty for 2D variables
# For single level: PRESSURE_LEVELS="850"
# For multiple levels (will be averaged): PRESSURE_LEVELS="850,500,300"
PRESSURE_LEVELS=""  # Pressure levels in hPa (empty for 2D variables)
# PRESSURE_LEVELS="500"  # Uncomment and set for 3D variables

# COF (Co-occurrence Feature) mask option
# COF_MASK=""  # Set to "--cof_mask" to extract COF masks instead of model variables
COF_MASK="--cof_mask"  # Uncomment to extract COF masks

# Extraction parameters
RADIUS="20.0"  # Extraction radius in degrees
LON_RES="0.25"  # Longitude resolution in degrees
LAT_RES="0.25"  # Latitude resolution in degrees

# Processing options
CHUNK_SIZE="1000"  # Chunk size for time dimension in zarr output
PROGRESS_FREQ="1000"  # How often to print progress

# Date filtering options (leave empty to process all tracks in the track file)
START_DATE=""  # e.g., "2019-08-01" or leave empty for no filtering
END_DATE=""    # e.g., "2020-08-31" or leave empty for no filtering

# Spatial filtering options (leave empty to use script defaults)
# Script defaults: min_lat = -90 + radius, max_lat = 90 - radius (avoids poles)
#                  min_lon/max_lon = None (no longitude filtering)
MIN_LAT=""   # e.g., "-80" or leave empty for script default
MAX_LAT=""   # e.g., "80" or leave empty for script default
MIN_LON=""   # e.g., "-180" or leave empty for no filtering
MAX_LON=""   # e.g., "180" or leave empty for no filtering

# Storm filtering (for testing - leave empty for production runs)
# STORM_IDS="11"  # Comma-separated storm IDs for testing
STORM_IDS=""  # Process all storms

# ===== SELECT VARIABLE FOR THIS ARRAY TASK =====
# Each array task processes one variable
CURRENT_VAR="${VARIABLES[$SLURM_ARRAY_TASK_ID]}"

echo "================================================"
echo "SLURM Array Job: $SLURM_ARRAY_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing variable: $CURRENT_VAR"
echo "Model: $CATALOG_MODEL (ERA5)"
echo "Grid type: Structured lat/lon"
echo "Track file: $TRACK_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "================================================"

# ===== BUILD COMMAND LINE ARGUMENTS =====
OPTIONAL_PARAMS=""

# Add date filtering if specified
if [ -n "$START_DATE" ] && [ -n "$END_DATE" ]; then
    OPTIONAL_PARAMS="$OPTIONAL_PARAMS --start_date $START_DATE --end_date $END_DATE"
    echo "Date range: $START_DATE to $END_DATE"
fi

# Add spatial filtering (only if specified)
if [ -n "$MIN_LAT" ]; then
    OPTIONAL_PARAMS="$OPTIONAL_PARAMS --min_lat $MIN_LAT"
fi
if [ -n "$MAX_LAT" ]; then
    OPTIONAL_PARAMS="$OPTIONAL_PARAMS --max_lat $MAX_LAT"
fi
if [ -n "$MIN_LON" ]; then
    OPTIONAL_PARAMS="$OPTIONAL_PARAMS --min_lon $MIN_LON"
fi
if [ -n "$MAX_LON" ]; then
    OPTIONAL_PARAMS="$OPTIONAL_PARAMS --max_lon $MAX_LON"
fi

# Add extraction parameters
OPTIONAL_PARAMS="$OPTIONAL_PARAMS --radius $RADIUS --lon_res $LON_RES --lat_res $LAT_RES"
OPTIONAL_PARAMS="$OPTIONAL_PARAMS --chunk_size $CHUNK_SIZE --progress_freq $PROGRESS_FREQ"

# Add pressure levels if specified (for 3D variables)
if [ -n "$PRESSURE_LEVELS" ]; then
    OPTIONAL_PARAMS="$OPTIONAL_PARAMS --pressure_levels $PRESSURE_LEVELS"
    echo "Pressure levels: $PRESSURE_LEVELS hPa"
fi

# Add storm ID filtering if specified (for testing)
if [ -n "$STORM_IDS" ]; then
    OPTIONAL_PARAMS="$OPTIONAL_PARAMS --storm_ids $STORM_IDS"
    echo "Testing mode: Processing only storm IDs: $STORM_IDS"
fi

# Add COF mask flag if specified
if [ -n "$COF_MASK" ]; then
    OPTIONAL_PARAMS="$OPTIONAL_PARAMS $COF_MASK"
    echo "Extracting COF masks"
fi

# IMPORTANT: Add --structured_mesh flag for ERA5 (lat/lon grid, not unstructured mesh)
OPTIONAL_PARAMS="$OPTIONAL_PARAMS --structured_mesh"

echo "================================================"
echo "Extraction Configuration:"
echo "  Model: $CATALOG_MODEL"
echo "  Grid type: Structured lat/lon (ERA5)"
echo "  Spatial bounds: [$MIN_LON, $MAX_LON] × [$MIN_LAT, $MAX_LAT]"
echo "  Extraction radius: ${RADIUS}°"
echo "  Grid resolution: ${LON_RES}° × ${LAT_RES}°"
echo "  Variable: $CURRENT_VAR"
if [ "$CURRENT_VAR" = "pr" ]; then
    echo "  NOTE: Using IMERG V7 precipitation for ERA5"
fi
if [ -n "$PRESSURE_LEVELS" ]; then
    echo "  Pressure levels: $PRESSURE_LEVELS hPa (for 3D variables)"
fi
echo "  Total variables in array: ${#VARIABLES[@]}"
echo "  Concurrent tasks limit: 3"
echo "================================================"

# ===== RUN EXTRACTION FOR SINGLE VARIABLE =====
# Process ONE variable per array task
# srun parameters:
#   -n 1: Run 1 task (single Python process for this variable)
#   -c 32: Allocate 32 CPU cores to the task (for Dask parallel operations)
#   --cpu_bind=cores: Bind threads to specific cores for better performance
# 
# Only the current variable (selected by SLURM_ARRAY_TASK_ID) is processed

# Build base command
CMD="srun -n 1 -c 32 --cpu_bind=cores python extract_etc_2d_vars.py \
  --catalog_url \"$CATALOG_URL\" \
  --catalog_model \"$CATALOG_MODEL\" \
  --catalog_params \"$CATALOG_PARAMS\" \
  --trackfile \"$TRACK_FILE\" \
  --output_dir \"$OUTPUT_DIR\" \
  --variables \"$CURRENT_VAR\" \
  $OPTIONAL_PARAMS"

# Add current_location only if not empty
if [ -n "$CURRENT_LOCATION" ]; then
    CMD="$CMD --current_location \"$CURRENT_LOCATION\""
fi

# Execute the command
eval $CMD

EXIT_CODE=$?

echo "================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: Variable $CURRENT_VAR completed at $(date)"
else
    echo "FAILED: Variable $CURRENT_VAR failed with exit code $EXIT_CODE at $(date)"
fi
echo "================================================"

exit $EXIT_CODE
