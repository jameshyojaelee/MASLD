#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=cohort-permutation
#SBATCH --array=1-200
#SBATCH --output=logs/cohort-permutation_%a.out
#SBATCH --error=logs/cohort-permutation_%a.err

# 10,000 total permutations: 200 array tasks x 50 perms each
set -eo pipefail
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq
set -u

export N_PERMS_PER_CHUNK=50
export B_TOTAL=10000

echo "=== cohort-permutation task ${SLURM_ARRAY_TASK_ID} ==="
echo "N_PERMS_PER_CHUNK=${N_PERMS_PER_CHUNK}  B_TOTAL=${B_TOTAL}  (total=10000)"
echo "Started: $(date)"

Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/14.5a_T3_cohort_perm.R

echo "Finished: $(date)"
