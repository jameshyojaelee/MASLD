#!/usr/bin/env bash
#SBATCH --job-name=agp_c1-endpoint2-gpu-recompute
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2_recompute_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2_recompute_%j.err
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PKG=${ROOT}/GWAS/finemapping/results/alphagenome_program/c1-endpoint2-gpu-20260914T203750Z
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

echo "host $(hostname) start $(date -u)"
"${PY}" "${ROOT}/scripts/analysis/alphagenome_program/c1e2_03_recompute.py" --outdir "${PKG}"
"${PY}" -m pip freeze > "${PKG}/pip_freeze.txt" 2>/dev/null || true
echo "done $(date -u)"
