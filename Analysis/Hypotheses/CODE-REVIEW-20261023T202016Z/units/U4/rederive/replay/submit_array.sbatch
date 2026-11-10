#!/bin/bash
#SBATCH --job-name=u4_replay_submit
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --dependency=afterok:22383429
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/logs/%j.out
# Reads the single-task timing from job 22383429, sizes memory (2x peak, 4-64 GB),
# sets the per-task time limit to 3x its wall time (>= 1 h), cuts the sample to
# 20 pairs if the task took > 3 h, and submits ONE array (%10).
set -eo pipefail
RP=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/replay
LOG=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/logs/22383429.out
K=$(awk '/^test task = /{print $4}' $LOG)
PEAK_KB=$(awk -F': ' '/Maximum resident set size/{print $2}' $LOG | tail -1)
WALL=$(awk -F': ' '/Elapsed \(wall clock\)/{print $2}' $LOG | tail -1)
SEC=$(echo "$WALL" | awk -F: '{s=0; for(i=1;i<=NF;i++) s=s*60+$i; printf "%d", s}')
TLIM_MIN=$(( (3 * SEC + 59) / 60 )); (( TLIM_MIN < 60 )) && TLIM_MIN=60   # 3x task 1, at least 1 h
MEM_GB=$(( (2 * PEAK_KB + 1048575) / 1048576 )); (( MEM_GB < 4 )) && MEM_GB=4; (( MEM_GB > 64 )) && MEM_GB=64
if (( SEC > 10800 )); then
  # > 3 h: keep 20 pairs, the first half of each stratum in task order
  IDX=$(awk -F'\t' 'NR>1{n[$2]++; t[$2,n[$2]]=$1; s[$2]=1} END{for(k in s) for(i=1;i<=int((n[k]+1)/2);i++) print t[k,i]}' $RP/pairs.tsv | sort -n | grep -vx "$K" | paste -sd,)
else
  IDX=$(awk -F'\t' 'NR>1{print $1}' $RP/pairs.tsv | grep -vx "$K" | paste -sd,)
fi
echo "timing task $K: peak ${PEAK_KB} kB, wall ${WALL} (${SEC} s) -> mem ${MEM_GB}G, time ${TLIM_MIN} min; array ${IDX}"
env -u SLURM_JOB_ID -u SLURM_JOBID sbatch --parsable --mem=${MEM_GB}G --time=${TLIM_MIN} --array=${IDX}%10 $RP/replay_array.sbatch | tee $RP/array_jobid.txt
