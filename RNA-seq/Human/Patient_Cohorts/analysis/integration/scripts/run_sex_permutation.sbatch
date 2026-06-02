#!/bin/bash
#SBATCH --job-name=sex_permute
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-50%20
#SBATCH --output=logs/sex_permute_%A_%a.out
#SBATCH --error=logs/sex_permute_%A_%a.err

mkdir -p logs

export PERM_ITER=$SLURM_ARRAY_TASK_ID

echo "Permutation iteration: ${PERM_ITER}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

Rscript 14.4a_sex_permutation_iter.R
