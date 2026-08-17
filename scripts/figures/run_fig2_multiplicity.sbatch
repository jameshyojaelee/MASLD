#!/bin/bash
#SBATCH --job-name=fig2mult
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig2mult_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig2mult_%j.log
set -euo pipefail
umask 027
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT="${ROOT}/GWAS/finemapping/results/fig2_multiplicity/$(date -u +%Y%m%dT%H%M%SZ)"
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "${ROOT}/logs" "${OUT}"
echo "[sbatch] out=${OUT} start=$(date -u)"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "${ROOT}/scripts/figures/fig2_multiplicity_correction.py" \
  --project-root "${ROOT}" --output-dir "${OUT}"
echo "[sbatch] DONE -> ${OUT}"
