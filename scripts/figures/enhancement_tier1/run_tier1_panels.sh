#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/tier1_panels_%j.out
#SBATCH --error=logs/tier1_panels_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/enhancement_tier1
mkdir -p logs
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
echo "[$(date)] generating Tier-1 enhancement panels on $(hostname)"
Rscript make_tier1_panels.R
echo "[$(date)] done"
