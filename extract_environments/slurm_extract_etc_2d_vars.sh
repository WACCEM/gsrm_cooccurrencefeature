#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -t 00:10:00
#SBATCH -J extract_scream
#SBATCH -A m1867
#SBATCH --mail-user=zhe.feng@pnnl.gov
#SBATCH --mail-type=FAIL,END

# ===== JOB STRUCTURE =====
# This script submits 1 SLURM job that runs 1 Python process
# All variables are processed sequentially in that single process
# Benefits:
#   - Track data loaded only once (saves time and memory)
#   - Efficient batched extraction (groups storms by timestamp)
#   - 32 cores available for internal parallelization (Dask operations)
# ===========================

# module load python
# module list
source activate /global/common/software/m1867/python/hackathon

# Set up paths and parameters
ROOT_DIR="/pscratch/sd/b/beharrop/kmscale_hackathon/hackathon_pre/screamv2_ne120_tracking"
TRACK_FILE="${ROOT_DIR}/screamv2_ne120_hp8.etc_stitched_nodes.txt"
OUTPUT_DIR="/pscratch/sd/w/wcmca1/hackathon/etc_data/tests/scream_ne120_inst/single_vars/"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# ===== PARAMETERS TO CUSTOMIZE =====
# Model and catalog settings
CATALOG_URL="https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml"
CURRENT_LOCATION="NERSC"
CATALOG_MODEL="scream_ne120_inst"  # Use scream_ne120_inst for instantaneous variables
CATALOG_PARAMS='{"zoom": 8}'

# Extraction parameters
RADIUS="10.0"  # Extraction radius in degrees
LON_RES="0.25"  # Longitude resolution in degrees
LAT_RES="0.25"  # Latitude resolution in degrees

# Processing options
CHUNK_SIZE="1000"  # Chunk size for time dimension in zarr output
PROGRESS_FREQ="1000"  # How often to print progress

# Set variables to extract (2D variables only, no pressure dimension)
# Common atmospheric variables:
VARIABLES=(
  "pr"      # precipitation
  "psl"     # sea level pressure
  "tas"     # surface air temperature
  "huss"    # surface specific humidity
  "clt"     # total cloud fraction
)

# Additional variables you might want:
# "ua850" "va850"  # 850 hPa winds (if available as 2D in catalog)
# "ua500" "va500"  # 500 hPa winds
# "zg500"          # 500 hPa geopotential height
# "rh850"          # 850 hPa relative humidity
# "uivt" "vivt"    # integrated vapor transport

# COF (Co-occurrence Feature) mask option
COF_MASK=""  # Set to "--cof_mask" to extract COF masks instead of model variables

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

echo "Starting ETC 2D variable extraction..."
echo "Processing ${#VARIABLES[@]} variable(s)"
echo "Track file: $TRACK_FILE"
echo "Output directory: $OUTPUT_DIR"

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
echo "  Variables: ${VARIABLES[@]}"
echo "================================================"

# ===== RUN EXTRACTION =====
# Process all variables in a SINGLE Python call (track data loaded once for all)
# srun parameters:
#   -n 1: Run 1 task (single Python process for all variables)
#   -c 32: Allocate 32 CPU cores to the task (for Dask parallel operations)
#   --cpu_bind=cores: Bind threads to specific cores for better performance
# 
# All variables in VARIABLES array are passed to Python at once.
# Python will:
#   1. Load track data once
#   2. Process each variable sequentially 
#   3. Each variable extraction can use the 32 cores internally for parallel operations
srun -n 1 -c 32 --cpu_bind=cores python extract_etc_2d_vars.py \
  --catalog_url "$CATALOG_URL" \
  --current_location "$CURRENT_LOCATION" \
  --catalog_model "$CATALOG_MODEL" \
  --catalog_params "$CATALOG_PARAMS" \
  --trackfile "$TRACK_FILE" \
  --output_dir "$OUTPUT_DIR" \
  --variables "${VARIABLES[@]}" \
  $OPTIONAL_PARAMS

echo "================================================"
echo "Extraction completed at $(date)"
echo "================================================"

# ===== OPTIONAL: COMBINE VARIABLES INTO SINGLE FILE =====
# Uncomment the following lines to automatically combine the extracted variables
# echo "Combining extracted variables into single zarr file..."
# python combine_etc_2d_vars.py \
#   --input_dir "$OUTPUT_DIR" \
#   --output_dir "$OUTPUT_DIR" \
#   --output_prefix "etc_2d_combined"

echo "All processing complete at $(date)"
