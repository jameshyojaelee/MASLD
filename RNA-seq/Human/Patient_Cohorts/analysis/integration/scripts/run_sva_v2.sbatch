#!/bin/bash
#SBATCH --job-name=SVA
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/SVA_v2_%j.out
#SBATCH --error=logs/SVA_v2_%j.err

# Pillar D — SVA v2: residuals-based sensitivity (kallisto rebuild).

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] Pillar D: SVA v2 (residuals) on $(hostname)"
Rscript dream_sva_residuals_v2.R
echo "[$(date)] Done SVA v2"
