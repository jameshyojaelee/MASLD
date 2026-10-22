#!/bin/bash
#SBATCH --job-name=6060
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H66-20261015T124530Z/slurm-6060-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H66-20261015T124530Z/slurm-6060-%j.err
set -euo pipefail
umask 027
H66=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H66-20261015T124530Z
OUT=${H66}/run-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}/tmp"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp MPLCONFIGDIR=${OUT}/tmp
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1   # ~/.local numexpr breaks pandas in rnaseq
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
${ENV}/bin/Rscript ${H66}/scripts/00_donor_map.R --out "${OUT}"
${ENV}/bin/python ${H66}/scripts/01_h66_observability.py --out "${OUT}" --seed 20260923
${ENV}/bin/python -m pip freeze > "${OUT}/pip_freeze.txt"
sha256sum ${H66}/scripts/* > "${OUT}/script_sha256.tsv"
