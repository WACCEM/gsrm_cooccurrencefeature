#!/bin/bash
# Bash script to run tasks sequentially on an interactive node
# Usage: 
#   1. Request an interactive node:
#      salloc -N 1 -C cpu -q interactive -t 04:00:00 -A m1867
#   2. Run this script:
#      bash /global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/run_interactive_IR_IMERG.sh

echo "=========================================="
echo "Starting sequential task execution for IR_IMERG"
echo "=========================================="
echo ""

# Activate Python environment
echo "Activating Python environment..."
source activate /global/common/software/m1867/python/pyflex-dev
echo ""

# Run each task sequentially
task_number=0
while IFS= read -r line; do
    task_number=$((task_number + 1))
    echo "=========================================="
    echo "Task ${task_number} of 4: $line"
    echo "=========================================="
    echo "Start time: $(date)"
    echo ""

    # Run the command
    $line

    # Check if command succeeded
    if [ $? -eq 0 ]; then
        echo ""
        echo "✓ Task ${task_number} completed successfully"
    else
        echo ""
        echo "✗ Task ${task_number} FAILED with exit code $?"
        echo "Stopping execution."
        exit 1
    fi

    echo "End time: $(date)"
    echo ""
done < /global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/tasks_all_IR_IMERG.txt

echo "=========================================="
echo "All 4 tasks completed successfully!"
echo "=========================================="
echo "Final end time: $(date)"
