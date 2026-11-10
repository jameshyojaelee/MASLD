#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=interactive
#SBATCH --partition=bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=500G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/worker/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/worker/%j.err
# SS S1 bigmem worker: s1_one.sh over the sample list, 8 libraries at a time (8 STAR threads each; measured STAR peak
# 39 GB, so about 312 GB at once; 500G is the bigmem partition minimum). Per-sample locks are shared with the nslab
# array. ORDER 1 starts at the top of the list, 2 at the bottom, 3 in the middle, so workers meet few locked samples.
set -uo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
F=$SS/scripts/s1_samples.txt; n=$(wc -l < $F)
case "${ORDER:-1}" in
  2) L=$(tac $F) ;;
  3) L=$( { tail -n +$((n/2+1)) $F; head -n $((n/2)) $F; } ) ;;
  *) L=$(cat $F) ;;
esac
echo "$L" | xargs -P 8 -I{} bash -c "bash $SS/scripts/s1_one.sh {} > $SS/logs/worker/${SLURM_JOB_ID}_{}.log 2>&1 || echo 'FAILED {}'"
echo "worker done"
