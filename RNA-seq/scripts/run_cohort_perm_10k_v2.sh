#!/bin/bash
#SBATCH --array=1-10000%50
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=96:00:00
#SBATCH --qos=nslab
#SBATCH --job-name=permutation
#SBATCH --output=RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/perm_%A_%a.out
#SBATCH --error=RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/perm_%A_%a.err

# Cohort permutation 10K v2 — weekend run
# 10,000 tasks × 1 perm each (serial, no BiocParallel)
# Fixes: 64G RAM (was 32G), 96h wall, idempotent skip
# Submit Friday: sbatch RNA-seq/scripts/run_cohort_perm_10k_v2.sh

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

export PERM_ITER=$SLURM_ARRAY_TASK_ID
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/dream_disease_permutation_iter.R
