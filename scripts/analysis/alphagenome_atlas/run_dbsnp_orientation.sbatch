#!/usr/bin/env bash
#SBATCH --job-name=ag_dbsnp_orient
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ag_dbsnp_orient_%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT="$1"
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd "$MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas"
$PY 68_dbsnp_orientation.py
echo "DBSNP_ORIENT_DONE"
