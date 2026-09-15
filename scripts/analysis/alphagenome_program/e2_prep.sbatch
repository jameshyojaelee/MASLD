#!/usr/bin/env bash
#SBATCH --job-name=agp_e2-transfer-prep
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e2-prep_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e2-prep_%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT="${ROOT}/GWAS/finemapping/results/alphagenome_program/e2-transfer-20260914T203752Z"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

"${PY}" "${ROOT}/scripts/analysis/alphagenome_program/e2_01_build_inputs.py" \
  --output "${OUT}/inputs"

printf 'agp e2 prep complete\n'
