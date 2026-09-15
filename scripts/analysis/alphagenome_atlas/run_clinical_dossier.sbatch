#!/usr/bin/env bash
#SBATCH --job-name=ag_p6d_dossier
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ag_p6d_dossier_%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT="$1"
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd "$MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas"
$PY 65_clinical_dossier.py
$PY 66_dossier_checks.py
echo "P6D_DOSSIER_DONE"
