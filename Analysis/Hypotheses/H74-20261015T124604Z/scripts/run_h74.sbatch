#!/bin/bash
#SBATCH --job-name=6140
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H74-20261015T124604Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H74-20261015T124604Z/logs/%j.err
set -euo pipefail
H74=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H74-20261015T124604Z
OUT=${H74}/run-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}/tmp"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp MPLCONFIGDIR=${OUT}/tmp
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
${ENV}/bin/python ${H74}/scripts/h74_nearest_gene.py --out "${OUT}" --seed 20260923
${ENV}/bin/python -m pip freeze > "${OUT}/pip_freeze.txt"
sha256sum ${H74}/scripts/* > "${OUT}/script_sha256.tsv"
