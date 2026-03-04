#!/bin/bash
#
# Wrapper script to run extract_mcs_cof_tracks.py for all sources.
#
# Usage: ./run_extract_mcs_cof_tracks_all.sh
#
# Author: Zhe Feng (zhe.feng@pnnl.gov)
# Date: 2026-03-04

# ---------------------------------------------------------------------------
# Source definitions: parallel arrays of source name and trackstats filename
# ---------------------------------------------------------------------------
SOURCES=(
    "scream"
    "icon_d3hp003"
    # "ifs_tco3999_rcbmf"
    "IMERGv7"
    "nicam_gl11"
    "casesm2_10km_nocumulus"
    "um_glm_n2560_RAL3p3"
)

# Corresponding MCS track stats filenames (same order as SOURCES)
TRACKSTATS_FILES=(
    "mcs_tracks_final_20190801.0000_20200901.0000.nc"
    "mcs_tracks_final_20200102.0000_20201231.2330.nc"
    # "mcs_tracks_final_20200101.0000_20210228.2330.nc"
    "mcs_tracks_final_20190801.0000_20200901.0000.nc"
    "mcs_tracks_final_20200301.0000_20210301.0000.nc"
    "mcs_tracks_final_20200301.0000_20210301.0000.nc"
    "mcs_tracks_final_20200201.0000_20210301.0000.nc"
)

# ---------------------------------------------------------------------------
# Base paths
# ---------------------------------------------------------------------------
MCS_ROOT="/pscratch/sd/w/wcmca1/hackathon/mcs"
COF_DIR="/pscratch/sd/w/wcmca1/hackathon/cof_masks"
OUTPUT_DIR="${COF_DIR}/stats"

# ---------------------------------------------------------------------------
# Script settings
# ---------------------------------------------------------------------------
COF_WINDOW=6
CHUNK_SIZE=100

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/extract_mcs_cof_tracks.py"

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
if [ ! -f "${PYTHON_SCRIPT}" ]; then
    echo "ERROR: Python script not found: ${PYTHON_SCRIPT}"
    exit 1
fi

if [ ${#SOURCES[@]} -ne ${#TRACKSTATS_FILES[@]} ]; then
    echo "ERROR: SOURCES and TRACKSTATS_FILES arrays must have the same length."
    exit 1
fi

# Activate conda environment
source activate /global/common/software/m1867/python/hackathon

# Create output directories
mkdir -p "${OUTPUT_DIR}"
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${LOG_DIR}"

MASTER_LOG="${LOG_DIR}/run_extract_mcs_cof_tracks_$(date +%Y%m%d_%H%M%S).log"

echo "========================================" | tee "${MASTER_LOG}"
echo "MCS COF Track Extraction - Batch Run"   | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"
echo "Start time    : $(date)"                 | tee -a "${MASTER_LOG}"
echo "Python script : ${PYTHON_SCRIPT}"        | tee -a "${MASTER_LOG}"
echo "MCS root      : ${MCS_ROOT}"             | tee -a "${MASTER_LOG}"
echo "COF dir       : ${COF_DIR}"              | tee -a "${MASTER_LOG}"
echo "Output dir    : ${OUTPUT_DIR}"           | tee -a "${MASTER_LOG}"
echo "COF window    : ${COF_WINDOW} h"         | tee -a "${MASTER_LOG}"
echo "Chunk size    : ${CHUNK_SIZE}"           | tee -a "${MASTER_LOG}"
echo "Num sources   : ${#SOURCES[@]}"          | tee -a "${MASTER_LOG}"
echo "Sources       : ${SOURCES[*]}"           | tee -a "${MASTER_LOG}"
echo "Master log    : ${MASTER_LOG}"           | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"
echo "" | tee -a "${MASTER_LOG}"

SUCCESS_COUNT=0
FAIL_COUNT=0

for i in "${!SOURCES[@]}"; do
    SOURCE="${SOURCES[$i]}"
    TRACKSTATS="${MCS_ROOT}/${SOURCE}/stats/${TRACKSTATS_FILES[$i]}"
    SOURCE_LOG="${LOG_DIR}/extract_mcs_cof_tracks_${SOURCE}_$(date +%Y%m%d_%H%M%S).log"

    echo "----------------------------------------" | tee -a "${MASTER_LOG}"
    echo "Source     : ${SOURCE}"                   | tee -a "${MASTER_LOG}"
    echo "Trackstats : ${TRACKSTATS}"               | tee -a "${MASTER_LOG}"
    echo "Start time : $(date)"                     | tee -a "${MASTER_LOG}"

    # Verify trackstats file exists before running
    if [ ! -f "${TRACKSTATS}" ]; then
        echo "✗ SKIPPED: trackstats file not found: ${TRACKSTATS}" | tee -a "${MASTER_LOG}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        echo "" | tee -a "${MASTER_LOG}"
        continue
    fi

    START_TIME=$(date +%s)

    python "${PYTHON_SCRIPT}" \
        --source "${SOURCE}" \
        --trackstats "${TRACKSTATS}" \
        --cof-dir "${COF_DIR}" \
        --output-dir "${OUTPUT_DIR}" \
        --cof-window "${COF_WINDOW}" \
        --chunk-size "${CHUNK_SIZE}" \
        2>&1 | tee "${SOURCE_LOG}"

    EXIT_CODE=${PIPESTATUS[0]}
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))

    if [ ${EXIT_CODE} -eq 0 ]; then
        echo "✓ SUCCESS: ${SOURCE} (${ELAPSED}s)" | tee -a "${MASTER_LOG}"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "✗ FAILED: ${SOURCE} (exit code: ${EXIT_CODE}, ${ELAPSED}s)" | tee -a "${MASTER_LOG}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi

    echo "End time : $(date)"        | tee -a "${MASTER_LOG}"
    echo "Log file : ${SOURCE_LOG}"  | tee -a "${MASTER_LOG}"
    echo "" | tee -a "${MASTER_LOG}"
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "========================================" | tee -a "${MASTER_LOG}"
echo "Batch Run Complete"                       | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"
echo "End time    : $(date)"                    | tee -a "${MASTER_LOG}"
echo "Total       : ${#SOURCES[@]}"             | tee -a "${MASTER_LOG}"
echo "Successful  : ${SUCCESS_COUNT}"           | tee -a "${MASTER_LOG}"
echo "Failed/Skip : ${FAIL_COUNT}"              | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"

echo "" | tee -a "${MASTER_LOG}"
echo "Output files:" | tee -a "${MASTER_LOG}"
ls -lh "${OUTPUT_DIR}"/*_mcs_cof_tracks_2d.nc    2>/dev/null | tee -a "${MASTER_LOG}"
ls -lh "${OUTPUT_DIR}"/*_mcs_cof_flags.parquet    2>/dev/null | tee -a "${MASTER_LOG}"
ls -lh "${OUTPUT_DIR}"/*_mcs_trackstats_cof.parquet 2>/dev/null | tee -a "${MASTER_LOG}"

if [ ${FAIL_COUNT} -gt 0 ]; then
    echo "" | tee -a "${MASTER_LOG}"
    echo "WARNING: ${FAIL_COUNT} source(s) failed or were skipped. Check logs in ${LOG_DIR}" | tee -a "${MASTER_LOG}"
    exit 1
fi

exit 0
