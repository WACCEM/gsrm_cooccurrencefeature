#!/bin/bash
#SBATCH --job-name=extremepcp
#SBATCH --account=m1867
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:20:00
#SBATCH --output=/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/log/extremepcp_%j.out
#SBATCH --error=/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/log/extremepcp_%j.err
#SBATCH --mail-user=zhe.feng@pnnl.gov
#SBATCH --mail-type=END

#
# SLURM batch script for extreme precipitation storm type attribution
# Submit with: sbatch slurm_stormtype_extreme_precip.sh [CATALOG_SOURCE] [PERCENTILES]
#
# Examples:
#   sbatch slurm_stormtype_extreme_precip.sh                      # Submits jobs for all default sources
#   sbatch slurm_stormtype_extreme_precip.sh scream_ne120 "P90 P95"
#   sbatch slurm_stormtype_extreme_precip.sh IR_IMERG "P90"
#

# Get arguments
CATALOG_SOURCE=${1}
PERCENTILES=${2:-"P90 P95"}

# Default sources to process if no source specified
DEFAULT_SOURCES=(
    "IR_IMERG"
    "scream_ne120"
    "icon_d3hp003"
    "nicam_gl11"
    "um_glm_n2560_RAL3p3"
    "casesm2_10km_nocumulus"
)

# Configuration
SCRIPT_DIR="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/scripts"
PYTHON_SCRIPT="${SCRIPT_DIR}/calc_stormtype_extreme_precip_spatial.py"
CONFIG_FILE="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources.yaml"
OUTPUT_DIR="/pscratch/sd/w/wcmca1/hackathon/extreme_precip"
LOG_DIR="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/log"

# Create logs directory
mkdir -p ${LOG_DIR}

# Check if this is a multi-source submission (no CATALOG_SOURCE provided)
if [ -z "${CATALOG_SOURCE}" ]; then
    echo "=========================================="
    echo "No CATALOG_SOURCE specified"
    echo "Submitting separate jobs for all default sources"
    echo "=========================================="
    echo ""
    
    for SOURCE in "${DEFAULT_SOURCES[@]}"; do
        echo "Submitting job for: ${SOURCE}"
        JOB_ID=$(sbatch --parsable \
            --job-name="extremepcp_${SOURCE}" \
            --output="${LOG_DIR}/extremepcp_${SOURCE}_%j.out" \
            --error="${LOG_DIR}/extremepcp_${SOURCE}_%j.err" \
            $0 ${SOURCE} "${PERCENTILES}")
        echo "  Job ID: ${JOB_ID}"
    done
    
    echo ""
    echo "=========================================="
    echo "✅ Submitted ${#DEFAULT_SOURCES[@]} jobs"
    echo "=========================================="
    exit 0
fi

# Single source execution below this point
# Load Python environment
module load python
source activate hackathon

# Set number of workers
N_WORKERS=${SLURM_CPUS_PER_TASK:-16}

# Set Dask configuration for optimal performance
export DASK_DISTRIBUTED__WORKER__MEMORY__TARGET=0.85
export DASK_DISTRIBUTED__WORKER__MEMORY__SPILL=0.90
export DASK_DISTRIBUTED__WORKER__MEMORY__PAUSE=0.95
export DASK_DISTRIBUTED__WORKER__MEMORY__TERMINATE=0.98

echo "=========================================="
echo "SLURM Job: Storm Type Spatial Attribution"
echo "=========================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Workers: ${N_WORKERS}"
echo "Catalog Source: ${CATALOG_SOURCE}"
echo "Percentiles: ${PERCENTILES}"
echo "=========================================="
echo ""

# Run the Python script
python ${PYTHON_SCRIPT} \
    --catalog_source ${CATALOG_SOURCE} \
    --config_file ${CONFIG_FILE} \
    --percentiles ${PERCENTILES} \
    --output_dir ${OUTPUT_DIR} \
    --n_workers ${N_WORKERS} \
    --compute_cloud_types

EXIT_CODE=$?

echo ""
echo "=========================================="
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✅ Job completed successfully!"
else
    echo "❌ Job failed with exit code: ${EXIT_CODE}"
fi
echo "=========================================="

exit ${EXIT_CODE}
