#!/bin/bash
#SBATCH --job-name=scanpy
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/fig4_spatial_candidates_%j.out
#SBATCH --error=logs/fig4_spatial_candidates_%j.err
# Render the full Fig 4 spatial candidate gallery (GSE192741 + Vu + CosMx).
# Output: figures/main/fig4_validation/spatial/_candidates/
# If cpu is busy, resubmit with: --qos=interactive (bypasses nslab 7TB cap; max 4 jobs)

set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
mkdir -p logs

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-local}"
eval "$(micromamba shell hook --shell bash)" 2>/dev/null || true
micromamba activate spatial

python scripts/figures/fig4_spatial_candidates.py
echo "[$(date)] finished rc=$?"
