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
# Submit with: sbatch slurm_extreme_precip_thresholds_1h.sh [CATALOG_SOURCE] [PERCENTILES] [AVAILABLE_MEMORY_GB] [ZOOM] [VERSION] [INPUT_ZARR]
#
# Examples:
#   sbatch slurm_extreme_precip_thresholds_1h.sh                       # Submits jobs for all 8 default sources (zoom 8)
#   sbatch slurm_extreme_precip_thresholds_1h.sh IMERG "95"
#   sbatch slurm_extreme_precip_thresholds_1h.sh GSMAP "95 99" 480
#   sbatch slurm_extreme_precip_thresholds_1h.sh scream_ne120 "90 95"
#   sbatch slurm_extreme_precip_thresholds_1h.sh GSMAP "95" 450 9 v1    # zoom 9, auto-subset to LONG_START/END_TIME
#   sbatch --time=08:00:00 slurm_extreme_precip_thresholds_1h.sh IMERG "95"   # overrides default 6h walltime for IMERG
#   sbatch --time=00:30:00 slurm_extreme_precip_thresholds_1h.sh casesm2_10km_nocumulus "95 99" 450 9 v1 \
#       /global/cfs/cdirs/wcm_shr/hk25/healpix/casesm2_10km_nocumulus/casesm2_10km_nocumulus_hourly_pr_and_rlut_hp9.zarr
#       # explicit INPUT_ZARR: overrides catalog loading entirely for any source (see --input_zarr in the python script)
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
# The only thing that legitimately differs by source is expected WALLTIME: IMERG and
# GSMAP take hours (multi-year); most model sources take minutes (~1 year); a few
# model sources (MEDIUM_SOURCES, currently um_glm_n2560_RAL3p3 and nicam_gl11 -- see
# their own comment below) need more than a fast source's default but nowhere near
# IMERG/GSMAP's budget, due to unusually fine native Zarr chunking. LONG_SOURCES/
# MEDIUM_SOURCES below exist purely to right-size --time per source in the
# multi-source submission loop; only LONG_SOURCES also gets a --start_time/--end_time
# subset -- the python command is otherwise identical for all sources.
#
# Zoom 9 sizing: measured from a real zoom-8, 24-year IMERG run
# (slurm/log/extremepcp1h_56434456.err: 3h21m total, ~2.75e7 values/s throughput,
# consistent between the all-period pass and the per-year loop). At zoom 9 over
# LONG_START_TIME..LONG_END_TIME (5 years, 4x the cells of zoom 8), the same throughput
# projects to roughly 2h50m total -- comfortably inside LONG_TIME's 6h default. If a
# per-year block ever approaches the node memory limit at zoom 9 (peak ~276GB at
# AVAILABLE_MEMORY_GB=450), lower AVAILABLE_MEMORY_GB (e.g. ~350) rather than raising
# --time; see determine_cell_chunk_size()'s docstring for why.

# Get arguments
CATALOG_SOURCE=${1}
PERCENTILES=${2:-"95 99"}
AVAILABLE_MEMORY_GB=${3:-450}
ZOOM=${4:-8}
VERSION=${5:-v1}
# Explicit local Zarr store override, for any source (see --input_zarr in the python
# script). Inherently source-specific, so this is only meaningful for a direct
# single-source submission below -- the bulk multi-source loop always passes exactly
# 5 positionals to each child job and so never sets this.
INPUT_ZARR=${6:-}

# Sources with multi-year records -- these are the only ones that need the long
# walltime, and the only ones subset to LONG_START_TIME/LONG_END_TIME below. The
# ~1-year model sources are always run over their full record. Update LONG_SOURCES if
# config_sources_1h.yaml gains another multi-year source.
LONG_SOURCES=("IMERG" "GSMAP")
LONG_TIME="06:00:00"

# Sources needing more than SHORT_TIME despite being ~1-year records. Measured after
# the lazy-chunking fix (_restore_lazy_chunking() in the python script) at zoom 9,
# 450GB budget: um_glm_n2560_RAL3p3 21m32s total, nicam_gl11 31m21s total -- the
# latter just past SHORT_TIME's 00:30:00, so it would be killed by a bulk submission
# without this tier. Both are here more for margin than measured need (SHORT_TIME
# would nearly cover um_glm already); add a source here if a similar close call shows
# up for it too.
MEDIUM_SOURCES=("um_glm_n2560_RAL3p3" "nicam_gl11")
MEDIUM_TIME="02:00:00"

SHORT_TIME="00:30:00"
LONG_START_TIME="2018-01-01T00"
LONG_END_TIME="2022-12-31T23"

# All sources to process if no CATALOG_SOURCE is specified
DEFAULT_SOURCES=(
    "IMERG"
    "GSMAP"
    "scream_ne120"
    "icon_d3hp003"
    "nicam_gl11"
    "um_glm_n2560_RAL3p3"
    "casesm2_10km_nocumulus"
    # "ifs_tco3999_rcbmf"
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

# True (exit 0) if $1 is in MEDIUM_SOURCES, false (exit 1) otherwise
is_medium_source() {
    local src=$1
    for s in "${MEDIUM_SOURCES[@]}"; do
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
        elif is_medium_source "${SOURCE}"; then
            TIME_LIMIT=${MEDIUM_TIME}
        else
            TIME_LIMIT=${SHORT_TIME}
        fi
        echo "Submitting job for: ${SOURCE} (time=${TIME_LIMIT}, zoom=${ZOOM})"
        JOB_ID=$(sbatch --parsable \
            --time="${TIME_LIMIT}" \
            --job-name="extremepcp1h_${SOURCE}" \
            --output="${LOG_DIR}/extremepcp1h_${SOURCE}_%j.out" \
            --error="${LOG_DIR}/extremepcp1h_${SOURCE}_%j.err" \
            $0 ${SOURCE} "${PERCENTILES}" ${AVAILABLE_MEMORY_GB} ${ZOOM} ${VERSION})
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
echo "Zoom level: ${ZOOM}"
echo "Version: ${VERSION}"
if [ -n "${INPUT_ZARR}" ]; then
    echo "Input Zarr (explicit override): ${INPUT_ZARR}"
fi

# Time subset applies only to the multi-year sources (IMERG/GSMAP); model sources
# always run over their full ~1-year record. Skipped entirely when INPUT_ZARR is set
# for a long source too -- that combination doesn't arise today (LONG_SOURCES is
# IMERG/GSMAP only, which don't need INPUT_ZARR), but if it ever does, --start_time/
# --end_time still apply on top of the explicit store, same as normal.
#
# Note: MEDIUM_SOURCES/MEDIUM_TIME above are not consulted anywhere below this point.
# A direct single-source submission (sbatch ... um_glm_n2560_RAL3p3 ...) always gets
# whatever --time was resolved at sbatch-parse-time -- the #SBATCH header's default
# (06:00:00, comfortably above MEDIUM_TIME's 2h) unless overridden with an explicit
# sbatch --time=... flag. MEDIUM_TIME only matters for the multi-source submission
# loop above, which explicitly passes --time="${TIME_LIMIT}" per child job.
TIME_ARGS=""
if is_long_source "${CATALOG_SOURCE}"; then
    TIME_ARGS="--start_time ${LONG_START_TIME} --end_time ${LONG_END_TIME}"
    echo "Time period: ${LONG_START_TIME} to ${LONG_END_TIME} (inclusive)"
else
    echo "Time period: full record (no subset applied for non-long sources)"
fi

INPUT_ZARR_ARGS=""
if [ -n "${INPUT_ZARR}" ]; then
    INPUT_ZARR_ARGS="--input_zarr ${INPUT_ZARR}"
fi
echo "=========================================="
echo ""

# Run the Python script.
# --cell_chunk_size is deliberately left unset so it auto-sizes per call (all-period vs.
# each year), aligned to the native Zarr chunk size and bounded by --available_memory_gb.
python ${PYTHON_SCRIPT} \
    --catalog_source ${CATALOG_SOURCE} \
    --config_file ${CONFIG_FILE} \
    --zoom ${ZOOM} \
    --percentiles ${PERCENTILES} \
    --min_precip_threshold 0.1 \
    --available_memory_gb ${AVAILABLE_MEMORY_GB} \
    ${TIME_ARGS} \
    ${INPUT_ZARR_ARGS} \
    --output_dir ${OUTPUT_DIR} \
    --version ${VERSION}

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
