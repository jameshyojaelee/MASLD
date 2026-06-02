#!/bin/bash
#SBATCH --job-name=tau-squared
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/M03b_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/M03b_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab

# M03b: Per-diet tau-squared classification
# Reviewer response: R2 #11, R3 P1-2

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

echo "=== M03b: Per-Diet Tau-Squared Classification ==="
echo "Start: $(date)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo ""

Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/scripts/M03b_per_diet_classification.R 2>&1

echo ""
echo "=== M03b complete: $(date) ==="
