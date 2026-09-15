#!/usr/bin/env bash
#SBATCH --job-name=fig5-promote
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=scripts/figures/logs/fig5_promotion_%j.out
#SBATCH --error=scripts/figures/logs/fig5_promotion_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
R_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "$PROJECT_ROOT"

"$PY_BIN" scripts/figures/promote_fig5_completion.py
"$R_BIN" Analysis/Multimodal_Program_Projection/scripts/06_validate_outputs.R
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage \
  figures/main/fig5_molecular_context/fig5_molecular_context.pdf
python3 scripts/manuscript/validate_resource_scope.py
printf 'PROMOTION_VALIDATED\t%s\n' \
  "figures/main/fig5_molecular_context/fig5_molecular_context.pdf"
