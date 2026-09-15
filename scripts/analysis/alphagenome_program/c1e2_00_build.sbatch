#!/usr/bin/env bash
#SBATCH --job-name=agp_c1-endpoint2-gpu-build
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2_build_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2_build_%j.err
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PKG=${ROOT}/GWAS/finemapping/results/alphagenome_program/c1-endpoint2-gpu-20260914T203750Z
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

echo "host $(hostname) start $(date -u)"
"${PY}" "${ROOT}/scripts/analysis/alphagenome_program/c1e2_00_build_variant_lists.py" \
  --outdir "${PKG}"
echo "done $(date -u)"
