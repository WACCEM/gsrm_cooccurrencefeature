#!/bin/bash
#
# Wrapper script to run ETC processing scripts for all sources:
# - combine_etc_2d_vars.py (combines individual zarr files)
# - calc_etc_spatial_stats.py
# - create_etc_composites.py
#
# Usage: ./run_calc_etc_spatial_stats_all.sh
#
# Author: Zhe Feng (zhe.feng@pnnl.gov)
# Date: 2026-02-17

# List of Python scripts to run (in order)
# Format: "relative_path/script_name.py" or just "script_name.py" for scripts in same directory
SCRIPTS=(
    "../extract_environments/combine_etc_2d_vars.py"
    "create_etc_composites.py"
    "calc_etc_spatial_stats.py"
)

# List of sources to process
SOURCES=(
    "era5"
    "scream"
    "icon_d3hp003"
    "um_glm_n2560_RAL3p3"
    "nicam_gl11"
    "casesm2_10km_nocumulus"
)

# Base paths
ZARR_PATH="/pscratch/sd/w/wcmca1/hackathon/etc_data"
OUTPUT_DIR="/pscratch/sd/w/wcmca1/hackathon/etc_data/stats"

# Script location (assume it's in the same directory as this wrapper)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if all Python scripts exist
for SCRIPT in "${SCRIPTS[@]}"; do
    # Resolve script path (supports relative paths)
    if [[ "${SCRIPT}" == /* ]]; then
        # Absolute path
        SCRIPT_PATH="${SCRIPT}"
    else
        # Relative path from SCRIPT_DIR
        SCRIPT_PATH="$(cd "${SCRIPT_DIR}" && cd "$(dirname "${SCRIPT}")" && pwd)/$(basename "${SCRIPT}")"
    fi
    
    if [ ! -f "${SCRIPT_PATH}" ]; then
        echo "ERROR: Python script not found: ${SCRIPT_PATH}"
        exit 1
    fi
done

# Activate conda environment (if needed)
source activate /global/common/software/m1867/python/hackathon

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Log file
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${LOG_DIR}"
MASTER_LOG="${LOG_DIR}/run_all_$(date +%Y%m%d_%H%M%S).log"

echo "========================================" | tee "${MASTER_LOG}"
echo "ETC Processing - Batch Processing" | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"
echo "Start time: $(date)" | tee -a "${MASTER_LOG}"
echo "Number of scripts: ${#SCRIPTS[@]}" | tee -a "${MASTER_LOG}"
echo "Scripts: ${SCRIPTS[*]}" | tee -a "${MASTER_LOG}"
echo "Number of sources: ${#SOURCES[@]}" | tee -a "${MASTER_LOG}"
echo "Sources: ${SOURCES[*]}" | tee -a "${MASTER_LOG}"
echo "Output directory: ${OUTPUT_DIR}" | tee -a "${MASTER_LOG}"
echo "Master log: ${MASTER_LOG}" | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"
echo "" | tee -a "${MASTER_LOG}"

# Counters
SUCCESS_COUNT=0
FAIL_COUNT=0
TOTAL_COUNT=$((${#SCRIPTS[@]} * ${#SOURCES[@]}))

# Process each script
for SCRIPT in "${SCRIPTS[@]}"; do
    SCRIPT_NAME=$(basename "${SCRIPT}" .py)
    # Resolve script path (supports relative paths)
    if [[ "${SCRIPT}" == /* ]]; then
        # Absolute path
        PYTHON_SCRIPT="${SCRIPT}"
    else
        # Relative path from SCRIPT_DIR
        PYTHON_SCRIPT="$(cd "${SCRIPT_DIR}" && cd "$(dirname "${SCRIPT}")" && pwd)/$(basename "${SCRIPT}")"
    fi
    
    echo "========================================" | tee -a "${MASTER_LOG}"
    echo "Running script: ${SCRIPT}" | tee -a "${MASTER_LOG}"
    echo "========================================" | tee -a "${MASTER_LOG}"
    echo "" | tee -a "${MASTER_LOG}"
    
    # Process each source for this script
    for SOURCE in "${SOURCES[@]}"; do
        echo "----------------------------------------" | tee -a "${MASTER_LOG}"
        echo "Script: ${SCRIPT} | Source: ${SOURCE}" | tee -a "${MASTER_LOG}"
        echo "Start time: $(date)" | tee -a "${MASTER_LOG}"
        
        # Individual log file for this script-source combination
        SOURCE_LOG="${LOG_DIR}/${SCRIPT_NAME}_${SOURCE}_$(date +%Y%m%d_%H%M%S).log"
        
        # Run the Python script with appropriate arguments
        START_TIME=$(date +%s)
        
        # Determine which arguments to use based on the script
        if [[ "${SCRIPT}" == *"combine_etc_2d_vars.py" ]]; then
            python "${PYTHON_SCRIPT}" \
                --source "${SOURCE}" \
                2>&1 | tee "${SOURCE_LOG}"
        elif [[ "${SCRIPT}" == *"calc_etc_spatial_stats.py" ]]; then
            python "${PYTHON_SCRIPT}" \
                --source "${SOURCE}" \
                --zarr-path "${ZARR_PATH}" \
                --output-dir "${OUTPUT_DIR}" \
                --basic-radius 10.0 \
                --mask-radii 10.0 15.0 \
                --pr-radius 10.0 \
                2>&1 | tee "${SOURCE_LOG}"
        elif [[ "${SCRIPT}" == *"create_etc_composites.py" ]]; then
            python "${PYTHON_SCRIPT}" \
                --source "${SOURCE}" \
                --zarr-path "${ZARR_PATH}" \
                --out-dir "${OUTPUT_DIR}/${SOURCE}" \
                2>&1 | tee "${SOURCE_LOG}"
        else
            # Generic fallback for other scripts
            python "${PYTHON_SCRIPT}" \
                --source "${SOURCE}" \
                --zarr-path "${ZARR_PATH}" \
                --output-dir "${OUTPUT_DIR}" \
                2>&1 | tee "${SOURCE_LOG}"
        fi
        
        # Check exit status
        EXIT_CODE=${PIPESTATUS[0]}
        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))
        
        if [ ${EXIT_CODE} -eq 0 ]; then
            echo "✓ SUCCESS: ${SCRIPT} | ${SOURCE} (${ELAPSED}s)" | tee -a "${MASTER_LOG}"
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        else
            echo "✗ FAILED: ${SCRIPT} | ${SOURCE} (exit code: ${EXIT_CODE})" | tee -a "${MASTER_LOG}"
            FAIL_COUNT=$((FAIL_COUNT + 1))
        fi
        
        echo "End time: $(date)" | tee -a "${MASTER_LOG}"
        echo "Log file: ${SOURCE_LOG}" | tee -a "${MASTER_LOG}"
        echo "" | tee -a "${MASTER_LOG}"
    done
    
    echo "" | tee -a "${MASTER_LOG}"
done

# Final summary
echo "========================================" | tee -a "${MASTER_LOG}"
echo "Batch Processing Complete" | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"
echo "End time: $(date)" | tee -a "${MASTER_LOG}"
echo "Total jobs (scripts × sources): ${TOTAL_COUNT}" | tee -a "${MASTER_LOG}"
echo "Successful: ${SUCCESS_COUNT}" | tee -a "${MASTER_LOG}"
echo "Failed: ${FAIL_COUNT}" | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"

# List output files
echo "" | tee -a "${MASTER_LOG}"
echo "Output files:" | tee -a "${MASTER_LOG}"
echo "Combined zarr files:" | tee -a "${MASTER_LOG}"
for SOURCE in "${SOURCES[@]}"; do
    if [ -d "${ZARR_PATH}/${SOURCE}" ]; then
        ls -d "${ZARR_PATH}/${SOURCE}"/etc_2d_combined_*.zarr 2>/dev/null | tee -a "${MASTER_LOG}"
    fi
done
echo "" | tee -a "${MASTER_LOG}"
echo "Statistics files:" | tee -a "${MASTER_LOG}"
ls -lh "${OUTPUT_DIR}"/etc_spatial_stats_*.nc 2>/dev/null | tee -a "${MASTER_LOG}"
echo "" | tee -a "${MASTER_LOG}"
echo "Composite files:" | tee -a "${MASTER_LOG}"
for SOURCE in "${SOURCES[@]}"; do
    if [ -d "${OUTPUT_DIR}/${SOURCE}" ]; then
        ls -lh "${OUTPUT_DIR}/${SOURCE}"/etc_composites_*.nc 2>/dev/null | tee -a "${MASTER_LOG}"
    fi
done

# Exit with error if any failed
if [ ${FAIL_COUNT} -gt 0 ]; then
    echo "" | tee -a "${MASTER_LOG}"
    echo "WARNING: ${FAIL_COUNT} source(s) failed. Check individual logs in ${LOG_DIR}" | tee -a "${MASTER_LOG}"
    exit 1
fi

exit 0
