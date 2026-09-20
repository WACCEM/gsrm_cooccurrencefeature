#!/bin/bash
# Run Analyses 1 and 2 of the COF pipeline for any set of sources on an interactive node (no queue wait, 4 h limit).
# It runs scripts/run_cof_pipeline.py in the background with nohup, so a lost connection does not stop it.
#
# Usage:
#   1. Request an interactive node (and log in to it if salloc does not put you there):
#        salloc -N 1 -C cpu -q interactive -t 04:00:00 -A m1867
#   2. Run this script with the options of run_cof_pipeline.py; --data-root is required (a test area, or the production
#      root /pscratch/sd/w/wcmca1/hackathon/ once you decide to promote):
#        bash slurm/run_interactive_cof_pipeline.sh --data-root /pscratch/sd/w/wcmca1/hackathon/tmp/round2 --dry-run
#        bash slurm/run_interactive_cof_pipeline.sh --data-root /pscratch/sd/w/wcmca1/hackathon/tmp/round2
#        bash slurm/run_interactive_cof_pipeline.sh --data-root DIR --sources scream icon --analysis 2
#        bash slurm/run_interactive_cof_pipeline.sh --data-root DIR --resume        # continue after a stop or a timeout
#      Add --foreground to keep it in the terminal (Ctrl-C stops the steps and leaves them resumable).
#   3. Watch it:
#        cat <data-root>/pipeline_logs/<run id>/status.json     # what is running, done, failed
#        tail -f <data-root>/pipeline_logs/launch_<time>.out    # the runner's events
#
# The expected wall time for all six sources is about 1.5 h (1.2-2.1 h); dry-run shows the schedule. Options and behaviour:
#   python scripts/run_cof_pipeline.py --help      and      docs/procedures/run_cof_pipeline.md

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${COF_PYTHON:-/global/common/software/m1867/python/hackathon/bin/python}"
RUNNER="${REPO}/scripts/run_cof_pipeline.py"

# The wrapper's own option, and the options that only print something (they run in the foreground)
FOREGROUND=0
ARGS=()
DATA_ROOT=""
PREV=""
for a in "$@"; do
    case "$a" in
        --foreground) FOREGROUND=1 ;;
        --dry-run|--preflight-only|--self-test|--help|-h) FOREGROUND=1; ARGS+=("$a") ;;
        *) ARGS+=("$a") ;;
    esac
    if [ "$PREV" = "--data-root" ]; then DATA_ROOT="$a"; fi
    PREV="$a"
done

if [ ! -x "$PYTHON" ]; then
    echo "Python not found: $PYTHON (set COF_PYTHON to the environment to use)"; exit 2
fi

if [ "$FOREGROUND" = "1" ]; then
    exec "$PYTHON" "$RUNNER" "${ARGS[@]}"
fi

if [ -z "$DATA_ROOT" ]; then
    echo "--data-root is required (a test area, or /pscratch/sd/w/wcmca1/hackathon/ for production)"; exit 2
fi

# Refuse a run that would fail at once (bad options, missing prerequisites, outputs that would be overwritten): the runner's
# dry run checks all of that in a second, before the background run is started.
if ! "$PYTHON" "$RUNNER" "${ARGS[@]}" --dry-run > /dev/null; then
    echo "The runner refused these options; showing why:"
    "$PYTHON" "$RUNNER" "${ARGS[@]}" --dry-run | tail -15
    exit 2
fi
# and the same for the input checks (a few seconds), so that a missing or incomplete input is reported here
if ! "$PYTHON" "$RUNNER" "${ARGS[@]}" --preflight-only > /tmp/cof_preflight_$$.txt 2>&1; then
    echo "The input checks failed:"; grep -i "problem\|failed" /tmp/cof_preflight_$$.txt | head -10
    rm -f /tmp/cof_preflight_$$.txt
    exit 2
fi
rm -f /tmp/cof_preflight_$$.txt

LAUNCH_DIR="${DATA_ROOT%/}/pipeline_logs"
mkdir -p "$LAUNCH_DIR"
LAUNCH_LOG="${LAUNCH_DIR}/launch_$(date +%Y%m%d_%H%M%S).out"
echo "Starting the pipeline in the background on $(hostname) ..."
nohup setsid "$PYTHON" "$RUNNER" "${ARGS[@]}" > "$LAUNCH_LOG" 2>&1 &
PID=$!
echo "  process id : ${PID}"
echo "  events     : tail -f ${LAUNCH_LOG}"
echo "  status     : cat ${LAUNCH_DIR}/*/status.json   (the newest run directory)"
echo "  stop       : kill -TERM ${PID}      (steps stop, the run can be continued with --resume)"
