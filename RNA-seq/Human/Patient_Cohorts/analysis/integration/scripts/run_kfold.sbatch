#!/bin/bash
#SBATCH --job-name=k-fold
#SBATCH --partition=cpu
#SBATCH --array=1-60%25
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --output=logs/k-fold_%A_%a.out
#SBATCH --error=logs/k-fold_%A_%a.err

# Pillar A — cohort-stratified k-fold CV (kallisto rebuild).
# Iters 1..10  use K=10  (each fold ~10% test);
# iters 11..60 use K=50 (each fold ~2%  test).

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
export KFOLD_ITER=${SLURM_ARRAY_TASK_ID}
echo "[$(date)] K-fold iter ${KFOLD_ITER} on $(hostname)"
Rscript dream_kfold_iter.R
echo "[$(date)] Done k-fold iter ${KFOLD_ITER}"
