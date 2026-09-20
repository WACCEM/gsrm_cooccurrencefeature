#!/bin/bash
#SBATCH -A m1657
#SBATCH -J cof_pipeline_src
#SBATCH -t 02:00:00
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --array=1-6
#SBATCH --signal=B:TERM@300
#SBATCH --output=/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/log/log_cof_pipeline_%A_%a.log
#SBATCH --mail-type=END
#SBATCH --mail-user=zhe.feng@pnnl.gov
#
# Analyses 1 and 2, ONE SOURCE PER NODE: array task n runs the whole chain (Steps 1-3, monthly, thresholds, attribution) for the n-th source
# on its own node. One queue wait for all six, but six nodes at once. Use it when the single-node job (slurm_run_cof_pipeline.sh) does
# not fit, for example when the node is shared or one source needs more workers than its share of the node.
#
# Usage (DATA_ROOT is required):
#   sbatch --export=ALL,DATA_ROOT=/pscratch/sd/w/wcmca1/hackathon/tmp/round2 slurm/slurm_run_cof_pipeline_array.sh
#   sbatch --array=1,3 --export=ALL,DATA_ROOT=DIR slurm/slurm_run_cof_pipeline_array.sh      # only IMERGv7 and icon_d3hp003
#   sbatch --export=ALL,DATA_ROOT=DIR,PIPELINE_ARGS="--resume" slurm/slurm_run_cof_pipeline_array.sh
# Order of the array index: 1 IMERGv7, 2 scream, 3 icon_d3hp003, 4 nicam_gl11, 5 um_glm_n2560_RAL3p3, 6 casesm2_10km_nocumulus.
# The state markers are per source, so the tasks share one data root without conflicts.

if [ -z "${DATA_ROOT:-}" ]; then
    echo "DATA_ROOT is not set (sbatch --export=ALL,DATA_ROOT=...)"; exit 2
fi
SOURCES_LIST=(IMERGv7 scream icon_d3hp003 nicam_gl11 um_glm_n2560_RAL3p3 casesm2_10km_nocumulus)
SRC=${SOURCES_LIST[$((SLURM_ARRAY_TASK_ID - 1))]}
date
echo "Array task ${SLURM_ARRAY_TASK_ID}: ${SRC}"
REPO=/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature
PYTHON=/global/common/software/m1867/python/hackathon/bin/python

# shellcheck disable=SC2086
"$PYTHON" "$REPO/scripts/run_cof_pipeline.py" --data-root "$DATA_ROOT" --sources "$SRC" --stagger-sec 10 ${PIPELINE_ARGS:-} &
RUNNER_PID=$!
trap 'kill -TERM $RUNNER_PID 2>/dev/null' TERM USR1
wait $RUNNER_PID
RC=$?
# a forwarded signal interrupts the first wait (status above 128); wait again for the runner to finish stopping its steps
if [ $RC -gt 128 ]; then wait $RUNNER_PID; RC=$?; fi
date
exit $RC
