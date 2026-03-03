    #!/bin/bash
#
# Bash script to calculate extreme precipitation percentiles for all data sources
# defined in config_sources_1h.yaml (1-hourly input data)
#
# Usage: bash run_all_extreme_precip_thresholds_1h.sh
#
# Author: Zhe Feng, zhe.feng@pnnl.gov
# Date: March 2026

# Configuration
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PYTHON_SCRIPT="${SCRIPT_DIR}/calc_extreme_precip_thresholds_1h.py"
CONFIG_FILE="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources_1h.yaml"

# List of all catalog sources from config_sources_1h.yaml
# Update this list if you add/remove sources in the config file
SOURCES=(
    "IMERG"
    "scream_ne120"
    "icon_d3hp003"
    "nicam_gl11"
    "um_glm_n2560_RAL3p3"
    "casesm2_10km_nocumulus"
    "ifs_tco3999_rcbmf"
)

# Optional: Customize these parameters as needed
# Uncomment and modify if you want different values than the defaults
# ZOOM=8
# PERCENTILES="90 95"
# TIME_DURATIONS="1h 6h 1D"
# METHOD="linear"
MIN_PRECIP_THRESHOLD=0.1  # Minimum precipitation threshold in mm/h (default: 0.1)
# OUTPUT_DIR="/pscratch/sd/w/wcmca1/hackathon/extreme_precip_1h/"
# VERSION="v1"

echo "========================================================================"
echo "Starting extreme precipitation percentile calculation for all sources"
echo "(1-hourly input data)"
echo "========================================================================"
echo "Python script: ${PYTHON_SCRIPT}"
echo "Config file: ${CONFIG_FILE}"
echo "Sources to process: ${SOURCES[@]}"
echo ""

# Loop through each source and run the Python script
for SOURCE in "${SOURCES[@]}"; do
    echo "========================================================================"
    echo "Processing: ${SOURCE}"
    echo "========================================================================"
    echo "Start time: $(date)"

    # Basic command with just the source (uses all defaults)
    CMD="python ${PYTHON_SCRIPT} --catalog_source ${SOURCE} --config_file ${CONFIG_FILE}"

    # Optional: Add custom parameters if defined above
    # if [ ! -z "${ZOOM}" ]; then
    #     CMD="${CMD} --zoom ${ZOOM}"
    # fi
    # if [ ! -z "${PERCENTILES}" ]; then
    #     CMD="${CMD} --percentiles ${PERCENTILES}"
    # fi
    # if [ ! -z "${TIME_DURATIONS}" ]; then
    #     CMD="${CMD} --time_durations ${TIME_DURATIONS}"
    # fi
    # if [ ! -z "${METHOD}" ]; then
    #     CMD="${CMD} --method ${METHOD}"
    # fi
    if [ ! -z "${MIN_PRECIP_THRESHOLD}" ]; then
        CMD="${CMD} --min_precip_threshold ${MIN_PRECIP_THRESHOLD}"
    fi
    # if [ ! -z "${OUTPUT_DIR}" ]; then
    #     CMD="${CMD} --output_dir ${OUTPUT_DIR}"
    # fi
    # if [ ! -z "${VERSION}" ]; then
    #     CMD="${CMD} --version ${VERSION}"
    # fi

    echo "Command: ${CMD}"
    echo ""
    
    # Execute the command
    eval ${CMD}

    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo ""
        echo "✓ Successfully completed: ${SOURCE}"
        echo "End time: $(date)"
    else
        echo ""
        echo "✗ Error processing: ${SOURCE}"
        echo "End time: $(date)"
        echo "Continuing with next source..."
    fi
    echo ""
done

echo "========================================================================"
echo "All sources processed!"
echo "Completion time: $(date)"
echo "========================================================================"
