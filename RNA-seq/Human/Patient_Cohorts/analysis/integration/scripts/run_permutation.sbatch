#!/bin/bash
#SBATCH --job-name=permutation-null
#SBATCH --partition=cpu
#SBATCH --array=1-1000%25
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --output=logs/permutation-null_%A_%a.out
#SBATCH --error=logs/permutation-null_%A_%a.err

# Pillar C — Permutation null (1000 label-shuffle dream fits).

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
export PERM_ITER=${SLURM_ARRAY_TASK_ID}
echo "[$(date)] Permutation iter ${PERM_ITER} on $(hostname)"
Rscript dream_disease_permutation_iter.R
echo "[$(date)] Done permutation iter ${PERM_ITER}"
