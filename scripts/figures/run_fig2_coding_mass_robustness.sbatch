#!/bin/bash
#SBATCH --job-name=fig2_robust
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig2_robust_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig2_robust_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
# new timestamped directory per run; never reuse 20260907T105813Z
OUT="${ROOT}/GWAS/finemapping/results/fig2_robustness/$(date -u +%Y%m%dT%H%M%SZ)"
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "${ROOT}/logs" "${OUT}"
echo "[sbatch] out=${OUT} start=$(date -u)"
P=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$P" "${ROOT}/scripts/figures/fig2_coding_mass_robustness.py" \
  --project-root "${ROOT}" --output-dir "${OUT}"
echo "FIG2_CODING_MASS_ROBUSTNESS_COMPLETE -> ${OUT}"
