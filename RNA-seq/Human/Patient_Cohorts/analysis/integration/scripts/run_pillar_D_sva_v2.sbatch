#!/bin/bash
#SBATCH --job-name=pillarD_sva2
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarD_sva2_%j.out
#SBATCH --error=logs/pillarD_sva2_%j.err

set -eo pipefail
mkdir -p logs

source ~/.bashrc 2>/dev/null || true
micromamba activate rnaseq

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] Starting SVA v2 (residuals-based) on $(hostname)"
Rscript dream_sva_residuals_v2.R
echo "[$(date)] Done SVA v2"
