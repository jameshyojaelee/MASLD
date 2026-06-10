#!/bin/bash
#SBATCH --job-name=lvqw
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/recompute_io_C2_%j.log
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
source ~/.bashrc 2>/dev/null || true
micromamba activate rnaseq
cd "$SLURM_SUBMIT_DIR"
set -eo pipefail
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/recompute_integration_only_cohort_support_C2.R
