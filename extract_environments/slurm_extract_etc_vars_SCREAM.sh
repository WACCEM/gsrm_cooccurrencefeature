#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
##SBATCH -q shared
##SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH -J scream
#SBATCH -A m1867
#SBATCH --array=0-9%10
#SBATCH --output=logs/extract_etc_SCREAM_%A_%a.log
#SBATCH --mail-user=zhe.feng@pnnl.gov
#SBATCH --mail-type=FAIL,END

# ===== JOB ARRAY STRUCTURE (SHARED QUEUE) =====
# This script uses SLURM job arrays on the SHARED queue
# --array=0-9%10 means:
#   - 10 tasks total (one per variable, indices 0-9)
#   - %10 limits to 10 simultaneous tasks (prevents overwhelming remote server)
# --mem=8G:
#   - Requests 8 GB memory per task (actual usage ~4-6 GB)
#   - Shared queue charges only for resources used
#   - More cost-efficient than full node (512 GB)
# Benefits:
#   - Cost: ~8 GB × 3 tasks = 24 GB vs 512 GB full node (~20x cheaper!)
#   - Parallel processing: 3x faster than sequential
#   - Fault tolerance: If one variable fails, others continue
#   - Resource isolation: Each variable gets dedicated resources
#   - Easy restart: Can resubmit only failed array indices
# 
# The concurrent setting can be overwriten at job submission:
#   sbatch --array=0-9%5 slurm_extract_etc_2d_vars.sh
# To restart failed tasks: 
#   sbatch --array=2,5 slurm_extract_etc_2d_vars.sh
# ================================

# module load python
# module list
source activate /global/common/software/m1867/python/hackathon

# Create logs directory if it doesn't exist
mkdir -p logs

# Set up paths and parameters
ROOT_DIR="/pscratch/sd/b/beharrop/kmscale_hackathon/hackathon_pre/screamv2_ne120_tracking"
TRACK_FILE="${ROOT_DIR}/screamv2_ne120_hp8.etc_stitched_nodes.txt"
OUTPUT_DIR="/pscratch/sd/w/wcmca1/hackathon/etc_data/scream/single_vars/"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# ===== PARAMETERS TO CUSTOMIZE =====
# Model and catalog settings
CATALOG_URL="https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml"
CURRENT_LOCATION="NERSC"
CATALOG_MODEL="scream_ne120_inst"  # Use scream_ne120_inst for 2D instantaneous variables
# CATALOG_MODEL="scream_ne120"  # Use scream_ne120 for 3D 3h average variables
# CATALOG_MODEL="scream2D_hrly"  # Use scream2D_hrly for 2D hourly variables
CATALOG_PARAMS='{"zoom": 8}'

# ===== VARIABLE CONFIGURATION =====
# 
# 2D VARIABLES (no pressure dimension):
#   pr, psl, ua850, va850, ua500, va500, rh850, uivt, vivt, zg500
#   mcs_ar_etc_overlap_mask, etc_mcs_ar_overlap_mask, ar_mcs_etc_overlap_mask
#
# 3D VARIABLES (require --pressure_levels):
#   ua, va, omega, hus (specific humidity)
#   Note: For 3D variables, specify pressure levels below
#
# Set variables to extract (this corresponds to SLURM array indices)
# E.g., for 4 variables, do: --array=0-3
VARIABLES=(
    "huss" "tas" "uas" "vas"
#   "pr" "psl" "ua850" "va850" "ua500" "va500" "rh850" "uivt" "vivt" "zg500"  # 2D variables
#   "ua" "va" "omega" "hus"  # 3D variables (need PRESSURE_LEVELS)
#   "mcs_ar_etc_overlap_mask"  # COF mask (2D, requires --cof_mask flag)
#   "ar_mcs_etc_overlap_mask"  # COF mask (2D, requires --cof_mask flag)
#   "etc_mcs_ar_overlap_mask"  # COF mask (2D, requires --cof_mask flag)
)

# 3D variable options (for pressure level data)
# Leave empty for 2D variables
# For single level: PRESSURE_LEVELS="850"
# For multiple levels (will be averaged): PRESSURE_LEVELS="850,500,300"
PRESSURE_LEVELS=""  # Pressure levels in hPa (empty for 2D variables)
# PRESSURE_LEVELS="500"  # Uncomment and set for 3D variables

# COF (Co-occurrence Feature) mask option
COF_MASK=""  # Set to "--cof_mask" to extract COF masks instead of model variables
# COF_MASK="--cof_mask"  # Set to "--cof_mask" to extract COF masks instead of model variables

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
# STORM_IDS="1014"  # Comma-separated storm IDs for testing
STORM_IDS=""  # Process all storms

# ===== SELECT VARIABLE FOR THIS ARRAY TASK =====
# Each array task processes one variable
CURRENT_VAR="${VARIABLES[$SLURM_ARRAY_TASK_ID]}"

echo "================================================"
echo "SLURM Array Job: $SLURM_ARRAY_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing variable: $CURRENT_VAR"
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

echo "================================================"
echo "Extraction Configuration:"
echo "  Model: $CATALOG_MODEL"
echo "  Spatial bounds: [$MIN_LON, $MAX_LON] × [$MIN_LAT, $MAX_LAT]"
echo "  Extraction radius: ${RADIUS}°"
echo "  Grid resolution: ${LON_RES}° × ${LAT_RES}°"
echo "  Variable: $CURRENT_VAR"
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
srun -n 1 -c 32 --cpu_bind=cores python extract_etc_2d_vars.py \
  --catalog_url "$CATALOG_URL" \
  --current_location "$CURRENT_LOCATION" \
  --catalog_model "$CATALOG_MODEL" \
  --catalog_params "$CATALOG_PARAMS" \
  --trackfile "$TRACK_FILE" \
  --output_dir "$OUTPUT_DIR" \
  --variables "$CURRENT_VAR" \
  $OPTIONAL_PARAMS

EXIT_CODE=$?

echo "================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: Variable $CURRENT_VAR completed at $(date)"
else
    echo "FAILED: Variable $CURRENT_VAR failed with exit code $EXIT_CODE at $(date)"
fi
echo "================================================"

exit $EXIT_CODE
