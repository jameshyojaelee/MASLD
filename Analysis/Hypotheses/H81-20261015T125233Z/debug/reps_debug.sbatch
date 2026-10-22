#!/bin/bash
#SBATCH --job-name=6213
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z/logs/%j.err
set -euo pipefail
umask 027
H81=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z
OUT=${H81}/debug/reps-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}/tmp" "${OUT}/replicates"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
sha256sum ${H81}/scripts/* > "${OUT}/script_sha256.tsv"
/usr/bin/time -v ${ENV}/bin/Rscript --vanilla ${H81}/scripts/01_h81_spikein.R --mode debug --out "${OUT}" --cores 3 --seed 20260923 --n-rep 3 --debug-unit-frac 0.5 --debug-genes-per-decile 700 --debug-other-genes 7000
