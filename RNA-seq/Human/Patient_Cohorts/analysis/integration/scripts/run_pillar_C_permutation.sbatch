#!/bin/bash
#SBATCH --job-name=pillarC_perm
#SBATCH --partition=cpu
#SBATCH --array=1-1000%30
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarC_perm_%A_%a.out
#SBATCH --error=logs/pillarC_perm_%A_%a.err

# Pillar C — within-cohort permutation null, B=1000.
# Shuffles group_binary within each dataset (preserves per-cohort case/control
# counts) and refits dream. Throttled to 30 concurrent (3 TB peak).

set -eo pipefail
mkdir -p logs

source ~/.bashrc 2>/dev/null || true
micromamba activate rnaseq

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
export PERM_ITER=${SLURM_ARRAY_TASK_ID}
echo "[$(date)] Starting permutation iter ${PERM_ITER} on $(hostname)"
Rscript dream_disease_permutation_iter.R
echo "[$(date)] Done permutation iter ${PERM_ITER}"
