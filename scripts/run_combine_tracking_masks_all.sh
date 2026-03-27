#!/bin/bash
#
# Run combine tracking masks for all sources
# Usage: bash run_combine_tracking_masks_all.sh [CATALOG_SOURCE]
#
# If no catalog source is provided, it processes all sources from the config file

# Configuration
SCRIPT_DIR="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/scripts"
PYTHON_SCRIPT="${SCRIPT_DIR}/combine_tracking_masks.py"
CONFIG_FILE="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources.yaml"

# Default: process all sources
if [ -z "$1" ]; then
    SOURCES=(
        "scream_ne120" 
        "icon_d3hp003" 
        "nicam_gl11" 
        "um_glm_n2560_RAL3p3" 
        "casesm2_10km_nocumulus"
    )
else
    SOURCES=("$1")
fi

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
        --source ${SOURCE} \
        --config ${CONFIG_FILE}
    
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