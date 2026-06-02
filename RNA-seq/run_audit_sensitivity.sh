#!/bin/bash
#SBATCH --job-name=liver_audit_sensitivity
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/audit_sensitivity_%j.log

set -euo pipefail

echo "=== Audit Sensitivity Analyses ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

Rscript 33_audit_sensitivity_analyses.R

echo "End: $(date)"
