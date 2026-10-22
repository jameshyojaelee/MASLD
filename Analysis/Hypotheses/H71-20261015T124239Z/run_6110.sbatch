#!/bin/bash
#SBATCH --job-name=6110
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H71-20261015T124239Z/slurm-6110-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H71-20261015T124239Z/slurm-6110-%j.err
set -euo pipefail
umask 027
H71=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H71-20261015T124239Z
OUT=${H71}/run-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}/tmp"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp MPLCONFIGDIR=${OUT}/tmp
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENVS=/gpfs/commons/home/jameslee/micromamba/envs
${ENVS}/spatial/bin/python ${H71}/scripts/01_xlsx_check.py > "${OUT}/xlsx_check.out"
${ENVS}/rnaseq/bin/Rscript ${H71}/scripts/00_control_expression.R --out "${OUT}"
${ENVS}/rnaseq/bin/python ${H71}/scripts/02_h71_requalify.py --out "${OUT}" --seed 20260923
${ENVS}/rnaseq/bin/python -m pip freeze > "${OUT}/pip_freeze.txt"
sha256sum ${H71}/scripts/* > "${OUT}/script_sha256.tsv"
