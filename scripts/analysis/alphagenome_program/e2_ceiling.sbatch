#!/usr/bin/env bash
#SBATCH --job-name=agp_e2-transfer-ceiling
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e2-ceiling_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e2-ceiling_%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT="${ROOT}/GWAS/finemapping/results/alphagenome_program/e2-transfer-20260914T203752Z"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

"${PY}" "${ROOT}/scripts/analysis/alphagenome_program/e2_04_label_ceiling.py" \
  --inputs "${OUT}/inputs" \
  --output "${OUT}" \
  --resamples 10000 \
  --bootstrap-seed 20260914

printf 'agp e2 measured-label ceiling complete\n'
