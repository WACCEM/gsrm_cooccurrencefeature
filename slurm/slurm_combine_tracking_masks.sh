#!/bin/bash
#SBATCH -A m1867
#SBATCH -J combine_tracking_masks
#SBATCH -t 00:05:00
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --output=/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/log/log_combine_tracking_masks_%A_%a.log
#SBATCH --mail-type=END
#SBATCH --mail-user=zhe.feng@pnnl.gov
#SBATCH --array=1-6

date
# Activate Python environment
source activate /global/common/software/m1867/python/pyflex-dev

# cd /global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/scripts/

# Takes a specified line ($SLURM_ARRAY_TASK_ID) from the task file
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm/tasks_combine_tracking_masks.txt)
echo $LINE
# Run the line as a command
$LINE

date
