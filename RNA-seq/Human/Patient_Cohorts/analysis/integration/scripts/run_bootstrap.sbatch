#!/bin/bash
#SBATCH --job-name=bootstrap
#SBATCH --partition=cpu
#SBATCH --array=1-1000%25
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --output=logs/bootstrap_%A_%a.out
#SBATCH --error=logs/bootstrap_%A_%a.err

# Pillar A — Bootstrap DEG stability (1000 iterations).

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
export BOOT_ITER=${SLURM_ARRAY_TASK_ID}
echo "[$(date)] Bootstrap iter ${BOOT_ITER} on $(hostname)"
Rscript dream_bootstrap_iter.R
echo "[$(date)] Done bootstrap iter ${BOOT_ITER}"
