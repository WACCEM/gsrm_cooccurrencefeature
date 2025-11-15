#!/bin/bash

# ===== INTERACTIVE RUN SCRIPT =====
# This script is for running ETC 2D extraction in an interactive terminal session
# 
# Usage:
#   1. Request an interactive node:
#      salloc -N 1 -C cpu -q interactive -t 04:00:00 -A m1867
#   2. Activate Python environment:
#      source activate /global/common/software/m1867/python/hackathon
#   3. Run this script:
#      bash run_extract_etc_vars_SCREAM.sh
#
# Benefits of interactive mode:
#   - Real-time output monitoring
#   - Easy debugging and testing
#   - Can interrupt and modify on the fly
# ===================================

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
# CATALOG_MODEL="scream_ne120_inst"  # Use scream_ne120_inst for 2D instantaneous variables
CATALOG_MODEL="scream_ne120"  # Use scream_ne120 for 3D 3h average variables
CATALOG_PARAMS='{"zoom": 8}'

# Extraction parameters
RADIUS="20.0"  # Extraction radius in degrees
LON_RES="0.25"  # Longitude resolution in degrees
LAT_RES="0.25"  # Latitude resolution in degrees

# Processing options
CHUNK_SIZE="1000"  # Chunk size for time dimension in zarr output
PROGRESS_FREQ="1000"  # How often to print progress

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
# Set variables to extract
VARIABLES=(
    "ua" "va" "omega"
    "hus"  # Example 3D variable - requires PRESSURE_LEVELS
#   "pr"  # Example 2D variable
#   "psl" "uivt" "vivt" "zg500"  # More 2D variables
#   "va" "omega" "hus"  # More 3D variables (need PRESSURE_LEVELS)
#   "mcs_ar_etc_overlap_mask"  # COF mask (2D, requires --cof_mask flag)
)

# 3D variable options (for pressure level data)
# Leave empty for 2D variables
# For single level: PRESSURE_LEVELS="850"
# For multiple levels (will be averaged among the layers): PRESSURE_LEVELS="800,750,700,600"
PRESSURE_LEVELS="850"  # Pressure levels in hPa
# PRESSURE_LEVELS=""  # Uncomment for 2D variables only

# COF (Co-occurrence Feature) mask option
COF_MASK=""
# COF_MASK="--cof_mask"  # Set to "--cof_mask" to extract COF masks instead of model variables

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
STORM_IDS="1014"  # Comma-separated storm IDs for testing
# STORM_IDS=""  # Uncomment to process all storms

echo "Starting ETC 2D variable extraction (Interactive Mode)..."
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
echo "  Variables: ${VARIABLES[@]}"
if [ -n "$PRESSURE_LEVELS" ]; then
    echo "  Pressure levels: $PRESSURE_LEVELS hPa (for 3D variables)"
fi
echo "================================================"

# ===== RUN EXTRACTION =====
# Interactive mode: Run directly with python (no srun)
# All variables are processed sequentially in a single Python process
# 
# Python will:
#   1. Load track data once
#   2. Process each variable sequentially 
#   3. Each variable extraction uses available cores for parallel operations
python extract_etc_2d_vars.py \
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
