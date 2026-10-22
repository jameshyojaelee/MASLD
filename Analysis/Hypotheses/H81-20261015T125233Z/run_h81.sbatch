#!/bin/bash
#SBATCH --job-name=6210
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z/logs/%j.err
set -euo pipefail
umask 027
H81=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z
OUT=${H81}/run-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}/tmp" "${OUT}/replicates"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
sha256sum ${H81}/scripts/* > "${OUT}/script_sha256.tsv"
# step 0-1: input hashes, guards G-A/G-B/G-C (stops on S1/S2; flags S3/S4)
${ENV}/bin/Rscript --vanilla ${H81}/scripts/01_h81_spikein.R --mode guard --out "${OUT}" --cores 11
# step 2: replicates 1..50 (each writes replicates/rep_NNN.tsv.gz plus a done marker once)
${ENV}/bin/Rscript --vanilla ${H81}/scripts/01_h81_spikein.R --mode replicates --out "${OUT}" --cores 16 --seed 20260923 --n-rep 50
# step 3: diagnostics first, then O1-O5 and the bootstrap
${ENV}/bin/Rscript --vanilla ${H81}/scripts/02_h81_summary.R --out "${OUT}" --seed 20261923
