#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/worker/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/worker/%j.err
# SS S1 sweep after the array and workers end: re-runs s1_one.sh on every library so any library missing a count
# (strand gate passed late) or genotypes (S0d list written late), or failed, is finished. STAR is redone only
# where the BAM or index is missing. Exits non-zero if any library is still incomplete.
set -uo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
if [ "${ORDER:-1}" == 2 ]; then L=$(tac $SS/scripts/s1_samples.txt); else L=$(cat $SS/scripts/s1_samples.txt); fi
echo "$L" | xargs -P 2 -I{} bash -c "bash $SS/scripts/s1_one.sh {} > $SS/logs/worker/${SLURM_JOB_ID}_{}.log 2>&1 || echo 'FAILED {}'"
echo "sweep pass done"
n_bad=0
A=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results/SteatoSITE/alignments/star
n_stop=0
while read -r s; do
  if [ ! -d "$A/$s" ] && compgen -G "$A/${s}_ambiguous_*" > /dev/null; then echo "stopped by mate-order probe: $s"; n_stop=$((n_stop+1)); continue; fi
  { { [ -e "$SS/done/$s" ] || [ -e "$SS/strand_fail_$s" ]; } && [ -e "$SS/geno_done/$s" ] && [ -e "$SS/out_s1d_ig/$s/ig_fragments.tsv" ]; } || { echo "incomplete $s"; n_bad=$((n_bad+1)); }
done < $SS/scripts/s1_samples.txt
echo "stopped by probe: $n_stop"
echo "incomplete libraries: $n_bad"; [ "$n_bad" -eq 0 ]
