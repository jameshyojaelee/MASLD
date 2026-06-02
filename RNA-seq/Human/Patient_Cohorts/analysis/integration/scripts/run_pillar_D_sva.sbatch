#!/bin/bash
#SBATCH --job-name=pillarD_sva
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarD_sva_%j.out
#SBATCH --error=logs/pillarD_sva_%j.err

# Pillar D — SVA-augmented dream sensitivity (housekeeping negative controls).
# Slightly higher memory budget for the svaseq SV factor estimation.

set -eo pipefail
mkdir -p logs

source ~/.bashrc 2>/dev/null || true
micromamba activate rnaseq

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] Starting SVA sensitivity on $(hostname)"
Rscript dream_sva_sensitivity.R
echo "[$(date)] Done SVA sensitivity"
