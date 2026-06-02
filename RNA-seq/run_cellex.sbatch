#!/bin/bash
#SBATCH --job-name=cellex
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/cellex_%j.out
#SBATCH --error=RNA-seq/logs/cellex_%j.err

# Phase 4B Step 1: CELLEX expression specificity
set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(micromamba shell hook -s bash)"
micromamba activate causal_inference

export PYTHONNOUSERSITE=1

echo "=== CELLEX Expression Specificity ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"

python RNA-seq/71_cellex_specificity.py

echo "=== Done: $(date) ==="
