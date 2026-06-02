#!/bin/bash
#SBATCH --job-name=pillarD_pvca
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarD_pvca_%j.out
#SBATCH --error=logs/pillarD_pvca_%j.err

# Pillar D — variancePartition on dream residuals.

set -eo pipefail
mkdir -p logs

source ~/.bashrc 2>/dev/null || true
micromamba activate rnaseq

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] Starting residual PVCA on $(hostname)"
Rscript dream_residual_pvca.R
echo "[$(date)] Done residual PVCA"
