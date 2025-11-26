#!/bin/bash
#
# Wrapper script to run calc_etc_spatial_stats.py for all sources
#
# Usage: ./run_calc_etc_spatial_stats_all.sh
#
# Author: Zhe Feng (zhe.feng@pnnl.gov)
# Date: 2025-11-25

# List of sources to process
SOURCES=(
    "era5"
    "scream"
    "icon_d3hp003"
    "um_glm_n2560_RAL3p3"
    "nicam_gl11"
)

# Base paths
ZARR_PATH="/pscratch/sd/w/wcmca1/hackathon/etc_data"
OUTPUT_DIR="/pscratch/sd/w/wcmca1/hackathon/etc_data/stats"

# Script location (assume it's in the same directory as this wrapper)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/calc_etc_spatial_stats.py"

# Check if Python script exists
if [ ! -f "${PYTHON_SCRIPT}" ]; then
    echo "ERROR: Python script not found: ${PYTHON_SCRIPT}"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Log file
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${LOG_DIR}"
MASTER_LOG="${LOG_DIR}/run_all_$(date +%Y%m%d_%H%M%S).log"

echo "========================================" | tee "${MASTER_LOG}"
echo "ETC Spatial Statistics - Batch Processing" | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"
echo "Start time: $(date)" | tee -a "${MASTER_LOG}"
echo "Number of sources: ${#SOURCES[@]}" | tee -a "${MASTER_LOG}"
echo "Sources: ${SOURCES[*]}" | tee -a "${MASTER_LOG}"
echo "Output directory: ${OUTPUT_DIR}" | tee -a "${MASTER_LOG}"
echo "Master log: ${MASTER_LOG}" | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"
echo "" | tee -a "${MASTER_LOG}"

# Counters
SUCCESS_COUNT=0
FAIL_COUNT=0
TOTAL_COUNT=${#SOURCES[@]}

# Process each source
for SOURCE in "${SOURCES[@]}"; do
    echo "----------------------------------------" | tee -a "${MASTER_LOG}"
    echo "Processing source: ${SOURCE}" | tee -a "${MASTER_LOG}"
    echo "Start time: $(date)" | tee -a "${MASTER_LOG}"
    
    # Individual log file for this source
    SOURCE_LOG="${LOG_DIR}/${SOURCE}_$(date +%Y%m%d_%H%M%S).log"
    
    # Run the Python script
    START_TIME=$(date +%s)
    
    python "${PYTHON_SCRIPT}" \
        --source "${SOURCE}" \
        --zarr-path "${ZARR_PATH}" \
        --output-dir "${OUTPUT_DIR}" \
        --basic-radius 10.0 \
        --mask-radii 10.0 15.0 \
        --pr-radius 10.0 \
        2>&1 | tee "${SOURCE_LOG}"
    
    # Check exit status
    EXIT_CODE=${PIPESTATUS[0]}
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    if [ ${EXIT_CODE} -eq 0 ]; then
        echo "✓ SUCCESS: ${SOURCE} (${ELAPSED}s)" | tee -a "${MASTER_LOG}"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "✗ FAILED: ${SOURCE} (exit code: ${EXIT_CODE})" | tee -a "${MASTER_LOG}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
    
    echo "End time: $(date)" | tee -a "${MASTER_LOG}"
    echo "Log file: ${SOURCE_LOG}" | tee -a "${MASTER_LOG}"
    echo "" | tee -a "${MASTER_LOG}"
done

# Final summary
echo "========================================" | tee -a "${MASTER_LOG}"
echo "Batch Processing Complete" | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"
echo "End time: $(date)" | tee -a "${MASTER_LOG}"
echo "Total sources: ${TOTAL_COUNT}" | tee -a "${MASTER_LOG}"
echo "Successful: ${SUCCESS_COUNT}" | tee -a "${MASTER_LOG}"
echo "Failed: ${FAIL_COUNT}" | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"

# List output files
echo "" | tee -a "${MASTER_LOG}"
echo "Output files:" | tee -a "${MASTER_LOG}"
ls -lh "${OUTPUT_DIR}"/etc_spatial_stats_*.nc 2>/dev/null | tee -a "${MASTER_LOG}"

# Exit with error if any failed
if [ ${FAIL_COUNT} -gt 0 ]; then
    echo "" | tee -a "${MASTER_LOG}"
    echo "WARNING: ${FAIL_COUNT} source(s) failed. Check individual logs in ${LOG_DIR}" | tee -a "${MASTER_LOG}"
    exit 1
fi

exit 0
