#!/bin/bash
#SBATCH --job-name=6273
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H87-20261015T130324Z/logs/slurm-6273-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H87-20261015T130324Z/logs/slurm-6273-%j.err
set -euo pipefail
umask 027
H87=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H87-20261015T130324Z
OUT=${H87}/debug-build-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}/tmp"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
${ENV}/bin/python ${H87}/scripts/02_build_coding_column.py --ensembl "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H87-20261015T130324Z/debug-fetch-20261015T130804Z/ensembl" --out "${OUT}" --seed 20260923 --debug-controls-only
