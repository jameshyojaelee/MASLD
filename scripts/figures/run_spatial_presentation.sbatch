#!/bin/bash
#SBATCH --job-name=spatial_panels
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/spatial_presentation_%j.out
#SBATCH --error=logs/spatial_presentation_%j.err

set -euo pipefail

echo "=== Spatial Presentation Panels ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/figures/logs

# Use rapids_singlecell env (has scanpy + anndata compatible versions)
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rapids_singlecell

# Ensure matplotlib can render without display
export MPLBACKEND=Agg

python scripts/figures/spatial_presentation_panels.py

echo ""
echo "=== DONE: $(date) ==="
