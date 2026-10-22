#!/bin/bash
#SBATCH --job-name=6260
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H86-20261015T125018Z/slurm-6260-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H86-20261015T125018Z/slurm-6260-%j.err
set -euo pipefail
umask 027
H86=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H86-20261015T125018Z
OUT=${H86}/run-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir "${OUT}"; mkdir "${OUT}/tmp"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
CAT=figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8
${ENV}/bin/Rscript ${H86}/scripts/01_resolve_symbols.R --out "${OUT}" --catalog "${CAT}" --seed 20260923
${ENV}/bin/python ${H86}/scripts/02_coverage.py --out "${OUT}" --catalog "${CAT}" --seed 20260923
${ENV}/bin/Rscript -e 'sessionInfo()' > "${OUT}/sessionInfo.txt"
${ENV}/bin/python -m pip freeze > "${OUT}/pip_freeze.txt"
sha256sum ${H86}/scripts/* > "${OUT}/script_sha256.tsv"
