#!/usr/bin/env bash
#SBATCH --job-name=fig5-final-check
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=scripts/figures/logs/fig5_final_validation_%j.out
#SBATCH --error=scripts/figures/logs/fig5_final_validation_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
R_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "$PROJECT_ROOT"

"$PY_BIN" scripts/figures/validate_fig5_canonical.py
"$R_BIN" Analysis/Multimodal_Program_Projection/scripts/06_validate_outputs.R
python3 scripts/manuscript/validate_resource_scope.py
