#!/bin/bash
#SBATCH --job-name=pillarB_xferr
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarB_xferr_%j.out
#SBATCH --error=logs/pillarB_xferr_%j.err

# Pillar B — cross-cohort disease prediction. Reuses existing dream LOO
# outputs in results/integration/loo_cv/ — does NOT refit dream — so this is
# light-weight (8 CPU, 60 GB) and fits well under the 7 TB cap.

set -eo pipefail
mkdir -p logs

source ~/.bashrc 2>/dev/null || true
micromamba activate rnaseq

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] Starting Pillar B transferability on $(hostname)"
Rscript pillar_B_transferability.R
echo "[$(date)] Done Pillar B"
