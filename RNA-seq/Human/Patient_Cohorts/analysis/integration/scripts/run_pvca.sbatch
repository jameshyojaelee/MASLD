#!/bin/bash
#SBATCH --job-name=PVCA
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/PVCA_%j.out
#SBATCH --error=logs/PVCA_%j.err

# Pillar D — variancePartition on dream residuals (kallisto rebuild).

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] Pillar D: PVCA on $(hostname)"
Rscript dream_residual_pvca.R
echo "[$(date)] Done PVCA"
