#!/bin/bash
#
# Run storm type spatial attribution for extreme precipitation
# Usage: bash run_stormtype_extreme_precip.sh [CATALOG_SOURCE]
#
# If no catalog source is provided, it processes all sources from the config file

# Configuration
SCRIPT_DIR="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/scripts"
PYTHON_SCRIPT="${SCRIPT_DIR}/calc_stormtype_extreme_precip_spatial.py"
CONFIG_FILE="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources.yaml"
OUTPUT_DIR="/pscratch/sd/w/wcmca1/hackathon/extreme_precip"

# Processing parameters
PERCENTILES="P95"
N_WORKERS=16  # Increase for better parallelization on NERSC

# Default: process all sources
if [ -z "$1" ]; then
    SOURCES=("scream_ne120" "IR_IMERG" "icon_d3hp003" "nicam_gl11" "um_glm_n2560_RAL3p3" "casesm2_10km_nocumulus")
else
    SOURCES=("$1")
fi

echo "=========================================="
echo "Storm Type Spatial Attribution Processing"
echo "=========================================="
echo "Catalog sources: ${SOURCES[@]}"
echo "Percentiles: ${PERCENTILES}"
echo "Workers: ${N_WORKERS}"
echo "Output: ${OUTPUT_DIR}"
echo "=========================================="
echo ""

# Process each source
for SOURCE in "${SOURCES[@]}"; do
    echo ""
    echo "=========================================="
    echo "Processing: ${SOURCE}"
    echo "=========================================="
    
    # Check if source exists in config
    if ! grep -q "^${SOURCE}:" "${CONFIG_FILE}"; then
        echo "⚠️  Warning: ${SOURCE} not found in config file"
        echo "Skipping..."
        continue
    fi
    
    # Run the Python script
    python ${PYTHON_SCRIPT} \
        --catalog_source ${SOURCE} \
        --config_file ${CONFIG_FILE} \
        --percentiles ${PERCENTILES} \
        --output_dir ${OUTPUT_DIR} \
        --n_workers ${N_WORKERS} \
        --compute_cloud_types
    
    if [ $? -eq 0 ]; then
        echo "✅ Successfully processed ${SOURCE}"
    else
        echo "❌ Error processing ${SOURCE}"
    fi
    
    echo ""
done

echo ""
echo "=========================================="
echo "All processing complete!"
echo "=========================================="
echo "Results saved in: ${OUTPUT_DIR}"
