#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=classifier-sensitivity
#SBATCH --output=logs/classifier-sensitivity_%j.out
#SBATCH --error=logs/classifier-sensitivity_%j.err

# Classifier sensitivity: N_BOOT=10000, N_PERM=10000 (hardcoded in script)
set -eo pipefail
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq
set -u

echo "=== classifier-sensitivity ==="
echo "N_BOOT=10000, N_PERM=10000 (set in 80_classifier_sensitivity.R)"
echo "Started: $(date)"

Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/80_classifier_sensitivity.R

echo "Finished: $(date)"
