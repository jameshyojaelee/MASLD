#!/bin/bash
#SBATCH --array=1-100%50
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=96:00:00
#SBATCH --qos=nslab
#SBATCH --job-name=bootstrap
#SBATCH --output=RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/boot_v2_%A_%a.out
#SBATCH --error=RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/boot_v2_%A_%a.err

# Bootstrap 10K v2 — weekend run
# 100 tasks × 100 reps = 10,000 iterations
# Fixes vs v1: 64G RAM (was 32G), 96h wall (was 48h), checkpoint every 10 reps
# Submit Friday evening: sbatch RNA-seq/scripts/run_bootstrap_10k_v2.sh

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

export BOOT_ITER=$SLURM_ARRAY_TASK_ID
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/dream_bootstrap_iter.R
