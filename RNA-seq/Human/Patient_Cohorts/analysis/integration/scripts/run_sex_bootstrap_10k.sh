#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=sex-bootstrap
#SBATCH --array=1-100
#SBATCH --output=logs/sex-bootstrap_%a.out
#SBATCH --error=logs/sex-bootstrap_%a.err

# 10,000 total bootstrap replicates: 100 array tasks x 100 reps each
set -eo pipefail
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq
set -u

export BOOT_ITER=${SLURM_ARRAY_TASK_ID}
export N_BOOT=100

echo "=== sex-bootstrap task ${SLURM_ARRAY_TASK_ID} ==="
echo "BOOT_ITER=${BOOT_ITER}  N_BOOT=${N_BOOT}  (total=10000)"
echo "Started: $(date)"

Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/14.4e_sex_bootstrap_T4.R

echo "Finished: $(date)"
