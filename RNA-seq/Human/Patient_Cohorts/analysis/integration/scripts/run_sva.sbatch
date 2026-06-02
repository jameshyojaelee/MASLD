#!/bin/bash
#SBATCH --job-name=SVA
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/SVA_%j.out
#SBATCH --error=logs/SVA_%j.err

# Pillar D — SVA-augmented dream sensitivity (kallisto rebuild).

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] Pillar D: SVA on $(hostname)"
Rscript dream_sva_sensitivity.R
echo "[$(date)] Done SVA"
