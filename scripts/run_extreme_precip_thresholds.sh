#!/bin/bash
#
# Flexible bash script to calculate extreme precipitation percentiles
# Can run one source or multiple sources
#
# Usage examples:
#   bash run_extreme_precip_thresholds.sh IR_IMERG
#   bash run_extreme_precip_thresholds.sh IR_IMERG scream_ne120 icon_d3hp003
#   bash run_extreme_precip_thresholds.sh --all
#
# Author: Zhe Feng, zhe.feng@pnnl.gov
# Date: November 2025

# Configuration
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PYTHON_SCRIPT="${SCRIPT_DIR}/calc_extreme_precip_thresholds.py"

# All available sources
ALL_SOURCES=(
    "IR_IMERG"
    "scream_ne120"
    "icon_d3hp003"
    "nicam_gl11"
    "um_glm_n2560_RAL3p3"
    "casesm2_10km_nocumulus"
    "ifs_tco3999_rcbmf"
)

# Check if no arguments provided
if [ $# -eq 0 ]; then
    echo "Error: No catalog source specified"
    echo ""
    echo "Usage: $0 SOURCE1 [SOURCE2 ...]"
    echo "   or: $0 --all"
    echo ""
    echo "Available sources:"
    for SOURCE in "${ALL_SOURCES[@]}"; do
        echo "  - ${SOURCE}"
    done
    echo ""
    echo "Examples:"
    echo "  $0 IR_IMERG"
    echo "  $0 IR_IMERG scream_ne120"
    echo "  $0 --all"
    exit 1
fi

# Determine which sources to process
if [ "$1" == "--all" ]; then
    SOURCES=("${ALL_SOURCES[@]}")
    echo "Processing all sources..."
else
    SOURCES=("$@")
    echo "Processing specified sources: ${SOURCES[@]}"
fi

echo "========================================================================"
echo "Starting extreme precipitation percentile calculation"
echo "========================================================================"
echo "Python script: ${PYTHON_SCRIPT}"
echo "Sources to process: ${SOURCES[@]}"
echo ""

# Loop through each source and run the Python script
SUCCESS_COUNT=0
FAIL_COUNT=0

for SOURCE in "${SOURCES[@]}"; do
    echo "========================================================================"
    echo "Processing: ${SOURCE}"
    echo "========================================================================"
    echo "Start time: $(date)"
    
    # Run with just the source (uses all defaults)
    CMD="python ${PYTHON_SCRIPT} --catalog_source ${SOURCE}"
    
    echo "Command: ${CMD}"
    echo ""
    
    # Execute the command
    eval ${CMD}
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo ""
        echo "✓ Successfully completed: ${SOURCE}"
        echo "End time: $(date)"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo ""
        echo "✗ Error processing: ${SOURCE}"
        echo "End time: $(date)"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        echo "Continuing with next source..."
    fi
    echo ""
done

echo "========================================================================"
echo "Processing complete!"
echo "========================================================================"
echo "Completion time: $(date)"
echo "Successfully processed: ${SUCCESS_COUNT} source(s)"
echo "Failed: ${FAIL_COUNT} source(s)"
echo "========================================================================"
