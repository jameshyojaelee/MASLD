#!/bin/bash
# Ascertainment control for the Figure 2 direct-trait coding-concentration claim.
# Streams ~29 GB of GWAS summary statistics to pull each credible set's lead
# variant effect size, so this must not run on the login node.
#SBATCH --job-name=ascertain
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ascertain_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ascertain_%j.log

set -euo pipefail
umask 027

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
STAMP=$(date -u +%Y%m%dT%H%M%SZ)
OUT="${ROOT}/GWAS/finemapping/results/ascertainment_control/${STAMP}"

export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

mkdir -p "${ROOT}/logs" "${OUT}"
echo "[sbatch] host=$(hostname) job=${SLURM_JOB_ID} out=${OUT} start=$(date -u)"
"${PY}" "${ROOT}/scripts/figures/fig2_ascertainment_control.py" \
  --project-root "${ROOT}" --output-dir "${OUT}"
echo "[sbatch] DONE $(date -u)  ->  ${OUT}"
