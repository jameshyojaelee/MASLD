#!/bin/bash
#SBATCH --job-name=6252
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H85-20261015T124643Z/logs/slurm-6252-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H85-20261015T124643Z/logs/slurm-6252-%j.err
# DEBUG run: earlier-set inputs of the outcome replaced by within-cohort permuted rows; NC-1 with 3 draws.
set -euo pipefail
umask 027
H85=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H85-20261015T124643Z
export H85_OUT=${H85}/debug-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${H85_OUT}"; mkdir -p "${H85_OUT}/tmp"
export TMPDIR=${H85_OUT}/tmp
export H85_DEBUG=1 H85_NC1_B=3
export PYTHONDONTWRITEBYTECODE=1
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export PATH=${ENV}/bin:${PATH}
sha256sum ${H85}/scripts/* > "${H85_OUT}/script_sha256.tsv"
${ENV}/bin/Rscript --vanilla ${H85}/scripts/program_axis_map_774_h85.R
