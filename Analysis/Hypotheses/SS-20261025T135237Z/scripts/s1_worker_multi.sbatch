#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=interactive
#SBATCH --partition=cpu,io
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/worker/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/worker/%j.err
# SS S1 multi-node worker: one job (one interactive-QOS slot) holding a 4 CPU / 48G piece on each of 10 nodes.
# Each node runs s1_one.sh serially over the library list, starting at a different offset; s1_one.sh skips
# finished and locked libraries, so pieces and other workers do not collide.
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
srun --ntasks=$SLURM_NNODES --ntasks-per-node=1 --cpus-per-task=4 --mem=48G --kill-on-bad-exit=0 bash -c '
  n=$(wc -l < '$SS'/scripts/s1_samples.txt); off=$(( (SLURM_PROCID * n) / SLURM_NTASKS ))
  { tail -n +$((off+1)) '$SS'/scripts/s1_samples.txt; head -n $off '$SS'/scripts/s1_samples.txt; } | while read -r S; do
    bash '$SS'/scripts/s1_one.sh "$S" > '$SS'/logs/worker/${SLURM_JOB_ID}_${SLURM_PROCID}_$S.log 2>&1 || echo "FAILED $S"
  done'
echo "multi worker done"
