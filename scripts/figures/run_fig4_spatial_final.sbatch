#!/bin/bash
#SBATCH --job-name=scanpy
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=logs/fig4_spatial_final_%j.out
#SBATCH --error=logs/fig4_spatial_final_%j.err
# Finalize the chosen panels into figures/main/fig4_validation/spatial/.

set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
mkdir -p logs
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-local}"
eval "$(micromamba shell hook --shell bash)" 2>/dev/null || true
micromamba activate spatial
python scripts/figures/fig4_spatial_final_panels.py
echo "[$(date)] finished rc=$?"
