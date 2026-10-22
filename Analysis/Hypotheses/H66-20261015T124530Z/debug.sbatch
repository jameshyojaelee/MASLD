#!/bin/bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
# DEBUG ONLY (not the prespec run). Usage: sbatch --job-name=NNNN --output=... debug.sbatch <mode> [extra args]
set -euo pipefail
umask 027
MODE=$1; shift
H66=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H66-20261015T124530Z
OUT=${H66}/debug-${MODE}-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}/tmp"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp MPLCONFIGDIR=${OUT}/tmp
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
${ENV}/bin/Rscript ${H66}/scripts/00_donor_map.R --out "${OUT}"
/usr/bin/time -v ${ENV}/bin/python ${H66}/scripts/01_h66_observability.py --out "${OUT}" --seed 20260923 --mode ${MODE} "$@"
