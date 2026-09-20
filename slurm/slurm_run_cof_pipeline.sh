#!/bin/bash
#SBATCH -A m1657
#SBATCH -J cof_pipeline
#SBATCH -t 03:00:00
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --signal=B:TERM@300
#SBATCH --output=/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/log/log_cof_pipeline_%j.log
#SBATCH --mail-type=END
#SBATCH --mail-user=zhe.feng@pnnl.gov
#
# Analyses 1 and 2 for all six sources (or a subset) in ONE job on ONE node: Steps 1-3, the monthly rain map, the extreme-precipitation
# thresholds and the attribution, in dependency order (scripts/run_cof_pipeline.py). One queue wait instead of one per step.
# Measured on a 503 GB node (all six sources, 2026-09-20): 95 min in total; the summed resident memory of all steps peaked at 266 GB while the six
# Step 1 runs overlapped. 3 h leaves about 2x margin.
#
# Usage (DATA_ROOT is required: a test area, or /pscratch/sd/w/wcmca1/hackathon/ for production):
#   sbatch --export=ALL,DATA_ROOT=/pscratch/sd/w/wcmca1/hackathon/tmp/round2 slurm/slurm_run_cof_pipeline.sh
#   sbatch --export=ALL,DATA_ROOT=DIR,SOURCES="scream icon",PIPELINE_ARGS="--analysis 2" slurm/slurm_run_cof_pipeline.sh
#   sbatch --export=ALL,DATA_ROOT=DIR,PIPELINE_ARGS="--resume" slurm/slurm_run_cof_pipeline.sh      # continue a stopped run
# SOURCES (default: all) and PIPELINE_ARGS (any other option of run_cof_pipeline.py) are optional.
# The job asks the runner to stop 5 minutes before the time limit (--signal), which leaves the steps resumable.

if [ -z "${DATA_ROOT:-}" ]; then
    echo "DATA_ROOT is not set (sbatch --export=ALL,DATA_ROOT=...)"; exit 2
fi
date
REPO=/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature
PYTHON=/global/common/software/m1867/python/hackathon/bin/python

SRC_ARGS=()
if [ -n "${SOURCES:-}" ]; then SRC_ARGS=(--sources ${SOURCES}); fi

# shellcheck disable=SC2086
"$PYTHON" "$REPO/scripts/run_cof_pipeline.py" --data-root "$DATA_ROOT" "${SRC_ARGS[@]}" ${PIPELINE_ARGS:-} &
RUNNER_PID=$!
# forward the time-limit warning and a cancel to the runner, which stops its steps cleanly
trap 'kill -TERM $RUNNER_PID 2>/dev/null' TERM USR1
wait $RUNNER_PID
RC=$?
# a forwarded signal interrupts the first wait (status above 128); wait again for the runner to finish stopping its steps
if [ $RC -gt 128 ]; then wait $RUNNER_PID; RC=$?; fi
date
exit $RC
