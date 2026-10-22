#!/bin/bash
#SBATCH --job-name=6270
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H87-20261015T130324Z/slurm-6270-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H87-20261015T130324Z/slurm-6270-%j.err
set -euo pipefail
umask 027
H87=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H87-20261015T130324Z
OUT=${H87}/run-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}/tmp" "${OUT}/ensembl"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
${ENV}/bin/python ${H87}/scripts/01_fetch_ensembl.py --out "${OUT}/ensembl"
${ENV}/bin/python ${H87}/scripts/02_build_coding_column.py --ensembl "${OUT}/ensembl" --out "${OUT}" --seed 20260923
${ENV}/bin/python -m pip freeze > "${OUT}/pip_freeze.txt"
sha256sum ${H87}/scripts/* > "${OUT}/script_sha256.tsv"
# Addition to the prespec header: copy the result to the run directory top level once.
test ! -e "${H87}/RESULT_H87.md" && cp "${OUT}/RESULT_H87.md" "${H87}/RESULT_H87.md"
echo "H87_DONE ${OUT}"
