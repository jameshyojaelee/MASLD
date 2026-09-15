#!/usr/bin/env bash
#SBATCH --job-name=ag_orient_check
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ag_orient_check_%j.out
set -euo pipefail
module load BCFtools/1.19-GCC-13.2.0
export PYTHONNOUSERSITE=1
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT="$1"
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd "$MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas"
$PY 68b_check_dbsnp_orientation.py
$PY -m pip freeze > "$AGA_OUT_ROOT/pip_freeze.txt"
echo "ORIENT_CHECK_DONE"
