#!/bin/bash
#SBATCH --job-name=chromVAR
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/chromvar_prop_%j.out
#SBATCH --error=logs/chromvar_prop_%j.err
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome
mkdir -p logs results/chromvar_v2
python scripts/03b_chromvar_per_donor_export.py
python scripts/09_chromvar_propagation.py
