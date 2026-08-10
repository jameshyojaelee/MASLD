#!/usr/bin/env bash
#SBATCH --job-name=spatial
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_%j.err

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PY_BIN="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
cd "$MASLD_PROJECT_ROOT"
"$PY_BIN" Analysis/Multimodal_Program_Projection/scripts/04_spatial_projection.py

