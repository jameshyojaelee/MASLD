#!/usr/bin/env bash
#SBATCH --job-name=ag_splice_detect
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ag_splice_detect_%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT="$1"
mkdir -p "$AGA_OUT_ROOT/tables"
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd "$MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas"
$PY -m unittest tests.test_splice_detection
$PY 35_splice_detection.py "$MASLD_PROJECT_ROOT/GWAS/finemapping/results/alphagenome_atlas/p3a-benchmarks-20260909T190214Z"
$PY -m pip freeze > "$AGA_OUT_ROOT/pip_freeze.txt"
echo "SPLICE_DETECTION_DONE"
