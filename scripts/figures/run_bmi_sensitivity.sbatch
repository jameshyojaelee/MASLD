#!/bin/bash
#SBATCH --job-name=bmi_sensitivity
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/bmi_sensitivity_%j.out
#SBATCH --error=logs/bmi_sensitivity_%j.err

set -euo pipefail

echo "=== BMI Confounding Sensitivity Analysis ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Started: $(date)"
echo ""

# Activate environment
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/figures/logs

Rscript scripts/figures/figS_bmi_sensitivity.R

echo ""
echo "Completed: $(date)"
