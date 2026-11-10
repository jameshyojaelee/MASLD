#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=interactive
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/worker/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/worker/%j.err
# SS S1 interactive-QOS worker (user allowed --qos=interactive, 2026-10-25): runs s1_one.sh over the sample list,
# 4 libraries at a time (8 threads, about 64G each). Order: forward for odd ORDER, reverse for even, so workers and
# the array start at different ends.
set -uo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
if [ "${ORDER:-1}" == 2 ]; then L=$(tac $SS/scripts/s1_samples.txt); else L=$(cat $SS/scripts/s1_samples.txt); fi
echo "$L" | xargs -P 4 -I{} bash -c "bash $SS/scripts/s1_one.sh {} > $SS/logs/worker/${SLURM_JOB_ID}_{}.log 2>&1 || echo 'FAILED {}'"
echo "worker done"
