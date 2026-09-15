#!/usr/bin/env bash
#SBATCH --job-name=agp_c2-prep
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c2-prep_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c2-prep_%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT="${ROOT}/GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

"${PY}" "${ROOT}/scripts/analysis/alphagenome_program/c2_01_build_variants.py" \
  --output "${OUT}/inputs"
