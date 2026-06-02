#!/bin/bash
#SBATCH --job-name=sex_subsample
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-250%20
#SBATCH --output=logs/sex_subsample_%A_%a.out
#SBATCH --error=logs/sex_subsample_%A_%a.err

mkdir -p logs

# Map SLURM_ARRAY_TASK_ID to (fraction, iteration):
# Tasks 1-50   -> frac=0.50, iter=1-50
# Tasks 51-100 -> frac=0.60, iter=1-50
# Tasks 101-150 -> frac=0.70, iter=1-50
# Tasks 151-200 -> frac=0.80, iter=1-50
# Tasks 201-250 -> frac=0.90, iter=1-50

FRACS=(0.50 0.60 0.70 0.80 0.90)
FRAC_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) / 50 ))
ITER=$(( (SLURM_ARRAY_TASK_ID - 1) % 50 + 1 ))
FRAC=${FRACS[$FRAC_IDX]}

echo "Task ${SLURM_ARRAY_TASK_ID}: FRAC=${FRAC}, ITER=${ITER}"

export SUBSAMPLE_FRAC=$FRAC
export SUBSAMPLE_ITER=$ITER

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

Rscript 14.4_sex_subsampling_iter.R
