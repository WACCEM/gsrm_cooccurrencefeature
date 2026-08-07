#!/bin/bash
#SBATCH --job-name=extremepcp1h
#SBATCH --account=m1867
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH --output=/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/log/extremepcp1h_%j.out
#SBATCH --error=/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/log/extremepcp1h_%j.err
#SBATCH --mail-user=zhe.feng@pnnl.gov
#SBATCH --mail-type=END,FAIL

#
# SLURM batch script for 1-hourly extreme precipitation percentiles
# (calc_extreme_precip_thresholds_1h.py): all-period percentiles, plus per-year
# percentiles and their interannual IQR for sources with multi-year records.
#
# Submit with: sbatch slurm_extreme_precip_thresholds_1h.sh [CATALOG_SOURCE] [PERCENTILES] [AVAILABLE_MEMORY_GB]
#
# Examples:
#   sbatch slurm_extreme_precip_thresholds_1h.sh                    # Submits jobs for all 8 default sources
#   sbatch slurm_extreme_precip_thresholds_1h.sh IMERG "95"
#   sbatch slurm_extreme_precip_thresholds_1h.sh GSMAP "95 99" 480
#   sbatch slurm_extreme_precip_thresholds_1h.sh scream_ne120 "90 95"
#   sbatch --time=08:00:00 slurm_extreme_precip_thresholds_1h.sh IMERG "95"   # overrides default 6h walltime for IMERG
#
# One command works for every source in config_sources_1h.yaml -- there is no
# IMERG/GSMAP-specific code path needed:
#   - calc_annual_precip_percentiles() requires >= 2 qualifying calendar years before
#     it computes anything annual/IQR-related; for the ~1-year model sources it logs
#     "Only 1 year(s) meet the minimum coverage..." and cleanly skips, leaving just the
#     all-period percentiles.
#   - determine_cell_chunk_size() sizes its processing block from the actual number of
#     time steps being processed, so a 1-year model's all-period pass is automatically
#     sized like a single year of IMERG/GSMAP's per-year loop (minutes, tens of GB).
# The only thing that legitimately differs by source is expected WALLTIME: IMERG (24 yr)
# and GSMAP (15 yr) take hours, everything else takes minutes. LONG_SOURCES below exists
# purely to request a shorter --time for the fast sources in the multi-source submission
# loop; the python command itself is identical for all sources.
#

# Get arguments
CATALOG_SOURCE=${1}
PERCENTILES=${2:-"95 99"}
AVAILABLE_MEMORY_GB=${3:-450}

# Sources with multi-year records -- these are the only ones that need the long walltime.
# Update this list if config_sources_1h.yaml gains another multi-year source.
LONG_SOURCES=("IMERG" "GSMAP")
LONG_TIME="06:00:00"
SHORT_TIME="00:30:00"

# All sources to process if no CATALOG_SOURCE is specified
DEFAULT_SOURCES=(
    "IMERG"
    "GSMAP"
    "scream_ne120"
    "icon_d3hp003"
    "nicam_gl11"
    "um_glm_n2560_RAL3p3"
    "casesm2_10km_nocumulus"
    "ifs_tco3999_rcbmf"
)

# Configuration
SCRIPT_DIR="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/scripts"
PYTHON_SCRIPT="${SCRIPT_DIR}/calc_extreme_precip_thresholds_1h.py"
CONFIG_FILE="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources_1h.yaml"
OUTPUT_DIR="/pscratch/sd/w/wcmca1/hackathon/extreme_precip_1h"
LOG_DIR="/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/log"

# Create logs directory
mkdir -p ${LOG_DIR}

# True (exit 0) if $1 is in LONG_SOURCES, false (exit 1) otherwise
is_long_source() {
    local src=$1
    for s in "${LONG_SOURCES[@]}"; do
        [ "$s" == "$src" ] && return 0
    done
    return 1
}

# Check if this is a multi-source submission (no CATALOG_SOURCE provided)
if [ -z "${CATALOG_SOURCE}" ]; then
    echo "=========================================="
    echo "No CATALOG_SOURCE specified"
    echo "Submitting separate jobs for all default sources"
    echo "=========================================="
    echo ""

    for SOURCE in "${DEFAULT_SOURCES[@]}"; do
        if is_long_source "${SOURCE}"; then
            TIME_LIMIT=${LONG_TIME}
        else
            TIME_LIMIT=${SHORT_TIME}
        fi
        echo "Submitting job for: ${SOURCE} (time=${TIME_LIMIT})"
        JOB_ID=$(sbatch --parsable \
            --time="${TIME_LIMIT}" \
            --job-name="extremepcp1h_${SOURCE}" \
            --output="${LOG_DIR}/extremepcp1h_${SOURCE}_%j.out" \
            --error="${LOG_DIR}/extremepcp1h_${SOURCE}_%j.err" \
            $0 ${SOURCE} "${PERCENTILES}" ${AVAILABLE_MEMORY_GB})
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

echo "=========================================="
echo "SLURM Job: 1-Hourly Extreme Precipitation Percentiles"
echo "=========================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Catalog Source: ${CATALOG_SOURCE}"
echo "Percentiles: ${PERCENTILES}"
echo "Available memory budget: ${AVAILABLE_MEMORY_GB} GB"
echo "=========================================="
echo ""

# Run the Python script.
# --cell_chunk_size is deliberately left unset so it auto-sizes per call (all-period vs.
# each year), aligned to the native Zarr chunk size and bounded by --available_memory_gb.
python ${PYTHON_SCRIPT} \
    --catalog_source ${CATALOG_SOURCE} \
    --config_file ${CONFIG_FILE} \
    --percentiles ${PERCENTILES} \
    --min_precip_threshold 0.1 \
    --available_memory_gb ${AVAILABLE_MEMORY_GB} \
    --output_dir ${OUTPUT_DIR} \
    --version v1

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
