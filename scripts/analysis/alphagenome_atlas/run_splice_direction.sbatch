#!/usr/bin/env bash
#SBATCH --job-name=ag_splice_direction
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ag_splice_direction_%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT="$1"
mkdir -p "$AGA_OUT_ROOT/tables"
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd "$MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas"
$PY -m unittest tests.test_splice_direction tests.test_splice_detection
$PY 36_splice_direction.py run "$MASLD_PROJECT_ROOT/GWAS/finemapping/results/alphagenome_atlas/p3a-benchmarks-20260909T190214Z" "${2:-0}"
$PY -m pip freeze > "$AGA_OUT_ROOT/pip_freeze.txt"
echo "SPLICE_DIRECTION_DONE"
