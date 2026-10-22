#!/bin/bash
#SBATCH --job-name=5700
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H57-20261015T112421Z/slurm-5700-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H57-20261015T112421Z/slurm-5700-%j.err
set -euo pipefail
umask 027
H57=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H57-20261015T112421Z
OUT=${H57}/run-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}/tmp"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp MPLCONFIGDIR=${OUT}/tmp
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
${ENV}/bin/Rscript ${H57}/scripts/00_control_expression.R --out "${OUT}"
${ENV}/bin/python ${H57}/scripts/01_h57_overlap.py --out "${OUT}" --seed 20260923
${ENV}/bin/python -m pip freeze > "${OUT}/pip_freeze.txt"
sha256sum ${H57}/scripts/* > "${OUT}/script_sha256.tsv"
