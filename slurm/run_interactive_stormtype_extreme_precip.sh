#!/bin/bash
# Bash script to run storm type extreme precip attribution sequentially on an interactive node
# Usage:
#   1. Request an interactive node:
#      salloc -N 1 -C cpu -q interactive -t 04:00:00 -A m1867
#   2. Run this script:
#      bash /global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/run_interactive_stormtype_extreme_precip.sh
#
# Optionally override the percentiles (default: "P90 P95"):
#      bash run_interactive_stormtype_extreme_precip.sh "P90 P95 P99"

# Percentiles to compute
PERCENTILES=${1:-"P95 P90"}

# Sources to process (same as slurm_stormtype_extreme_precip.sh DEFAULT_SOURCES)
SOURCES=(
    "IR_IMERG"
    "scream_ne120"
    "icon_d3hp003"
    "nicam_gl11"
    "um_glm_n2560_RAL3p3"
    "casesm2_10km_nocumulus"
)
N_SOURCES=${#SOURCES[@]}

# Configuration (matches slurm_stormtype_extreme_precip.sh)
SCRIPT_DIR="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/scripts"
PYTHON_SCRIPT="${SCRIPT_DIR}/calc_stormtype_extreme_precip_spatial.py"
CONFIG_FILE="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources.yaml"
OUTPUT_DIR="/pscratch/sd/w/wcmca1/hackathon/extreme_precip"
LOG_DIR="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/log"

# Create logs directory
mkdir -p ${LOG_DIR}

echo "=========================================="
echo "Starting sequential task execution for stormtype_extreme_precip"
echo "=========================================="
echo "Percentiles: ${PERCENTILES}"
echo "Sources (${N_SOURCES}): ${SOURCES[*]}"
echo "=========================================="
echo ""

# Activate Python environment
echo "Activating Python environment..."
module load python
source activate hackathon
echo ""

# Set number of workers
N_WORKERS=${SLURM_CPUS_PER_TASK:-16}

# Set Dask configuration for optimal performance
export DASK_DISTRIBUTED__WORKER__MEMORY__TARGET=0.85
export DASK_DISTRIBUTED__WORKER__MEMORY__SPILL=0.90
export DASK_DISTRIBUTED__WORKER__MEMORY__PAUSE=0.95
export DASK_DISTRIBUTED__WORKER__MEMORY__TERMINATE=0.98

# Run each source sequentially
task_number=0
for SOURCE in "${SOURCES[@]}"; do
    task_number=$((task_number + 1))
    echo "=========================================="
    echo "Task ${task_number} of ${N_SOURCES}: ${SOURCE}"
    echo "=========================================="
    echo "Start time: $(date)"
    echo ""

    # Run the command
    python ${PYTHON_SCRIPT} \
        --catalog_source ${SOURCE} \
        --config_file ${CONFIG_FILE} \
        --percentiles ${PERCENTILES} \
        --output_dir ${OUTPUT_DIR} \
        --n_workers ${N_WORKERS} \
        --compute_cloud_types

    EXIT_CODE=$?

    # Check if command succeeded
    if [ ${EXIT_CODE} -eq 0 ]; then
        echo ""
        echo "✓ Task ${task_number} (${SOURCE}) completed successfully"
    else
        echo ""
        echo "✗ Task ${task_number} (${SOURCE}) FAILED with exit code ${EXIT_CODE}"
        echo "Stopping execution."
        exit 1
    fi

    echo "End time: $(date)"
    echo ""
done

echo "=========================================="
echo "All ${N_SOURCES} tasks completed successfully!"
echo "=========================================="
echo "Final end time: $(date)"
