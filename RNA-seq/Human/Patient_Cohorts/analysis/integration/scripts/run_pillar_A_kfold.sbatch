#!/bin/bash
#SBATCH --job-name=pillarA_kfold
#SBATCH --partition=cpu
#SBATCH --array=1-60%30
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarA_kfold_%A_%a.out
#SBATCH --error=logs/pillarA_kfold_%A_%a.err

# Pillar A — cohort-stratified k-fold CV.
# Iters 1..10  use K=10  (each fold ~10% test);
# iters 11..60 use K=50 (each fold ~2%  test).

set -eo pipefail
mkdir -p logs

source ~/.bashrc 2>/dev/null || true
micromamba activate rnaseq

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
export KFOLD_ITER=${SLURM_ARRAY_TASK_ID}
echo "[$(date)] Starting kfold iter ${KFOLD_ITER} on $(hostname)"
Rscript dream_kfold_iter.R
echo "[$(date)] Done kfold iter ${KFOLD_ITER}"
