#!/bin/bash
#SBATCH --job-name=6221
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H82-20261015T125052Z/logs/slurm-6221-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H82-20261015T125052Z/logs/slurm-6221-%j.err
# H82 debug job (synth); never a result.
set -euo pipefail
umask 027
H82=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H82-20261015T125052Z
OUT=${H82}/debug/synth-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}/tmp"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
${ENV}/bin/Rscript ${H82}/scripts/01_h82_refits.R --out "${OUT}" --seed 20260923 --workers 16 --B 4 --synthetic
${ENV}/bin/Rscript ${H82}/scripts/02_h82_fig_s3_panel.R --out "${OUT}"
sha256sum ${H82}/scripts/* > "${OUT}/script_sha256.tsv"
