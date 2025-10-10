#!/bin/sh
#SBATCH -A m1867
#SBATCH -J casesm2_mcsmonthly
#SBATCH -N 2 -c 128
#SBATCH -q regular
#SBATCH -t 00:15:00
#SBATCH -C cpu
#SBATCH --mail-user=zhe.feng@pnnl.gov
#SBATCH --mail-type=END
#SBATCH --output=log_mcs_monthly_casesm2_10km_cumulus.log

date

cd /global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/slurm
export THREADS=128

runcommands.sh tasklist_mcs_monthly_casesm2_10km_cumulus.txt

date