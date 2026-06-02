#!/bin/bash
#SBATCH --job-name=p1_loo
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p1_loo_%x_%j.out
#SBATCH --error=logs/p1_loo_%x_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "Host: $(hostname)"
echo "Started: $(date) — HELD_OUT=${HELD_OUT}"

Rscript dream_loo_cv.R

echo "Finished: $(date)"
