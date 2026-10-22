#!/bin/bash
#SBATCH --job-name=6030
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H63-20261015T125909Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H63-20261015T125909Z/logs/%j.err
set -euo pipefail
H63=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H63-20261015T125909Z
OUT=${H63}/run-$(date -u +%Y%m%dT%H%M%SZ)
mkdir "${OUT}"; mkdir "${OUT}/tmp"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp MPLCONFIGDIR=${OUT}/tmp
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
sha256sum ${H63}/PRESPEC_H63.md ${H63}/scripts/* > "${OUT}/script_sha256.tsv"
${ENV}/bin/python ${H63}/scripts/run_h63.py --out "${OUT}" --seed 20260923
${ENV}/bin/python -m pip freeze > "${OUT}/pip_freeze.txt"
