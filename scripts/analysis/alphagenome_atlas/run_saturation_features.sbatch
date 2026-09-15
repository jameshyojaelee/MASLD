#!/usr/bin/env bash
#SBATCH --job-name=ag_sat_features
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ag_sat_features_%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT="$1"
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd "$MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas"
$PY 43_saturation_features.py "$MASLD_PROJECT_ROOT/GWAS/finemapping/results/alphagenome_atlas/p4-saturation-20260909T191917Z" "${2:-8}" "${3:-0}"
$PY -m pip freeze > "$AGA_OUT_ROOT/pip_freeze.txt"
echo "SAT_FEATURES_DONE"
