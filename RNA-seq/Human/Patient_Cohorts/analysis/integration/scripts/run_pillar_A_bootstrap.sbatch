#!/bin/bash
#SBATCH --job-name=pillarA_boot
#SBATCH --partition=cpu
#SBATCH --array=1-1000%30
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarA_boot_%A_%a.out
#SBATCH --error=logs/pillarA_boot_%A_%a.err

# Pillar A — cohort-stratified subsampling-without-replacement, B=1000.
# Throttled to 30 concurrent (3 TB peak) to share the 7 TB nslab cap with
# the parallel CPSS / kfold / permutation arrays.

set -eo pipefail
mkdir -p logs

source ~/.bashrc 2>/dev/null || true
micromamba activate rnaseq

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
export BOOT_ITER=${SLURM_ARRAY_TASK_ID}
echo "[$(date)] Starting bootstrap iter ${BOOT_ITER} on $(hostname)"
Rscript dream_bootstrap_iter.R
echo "[$(date)] Done bootstrap iter ${BOOT_ITER}"
