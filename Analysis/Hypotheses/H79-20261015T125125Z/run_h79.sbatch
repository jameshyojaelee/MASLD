#!/bin/bash
#SBATCH --job-name=6190
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H79-20261015T125125Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H79-20261015T125125Z/logs/%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
H79=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H79-20261015T125125Z
OUT=${H79}/results/run-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}/tmp"
export TMPDIR=${OUT}/tmp H79_OUT=${OUT} H79_CORES=16 H79_SEED=20260923
cp ${H79}/scripts/h79_ancestry_refit.R ${OUT}/h79_ancestry_refit.as_run.R
sha256sum ${H79}/scripts/* > ${OUT}/script_sha256.tsv
echo "job ${SLURM_JOB_ID} host $(hostname) start $(date -u)"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript ${H79}/scripts/h79_ancestry_refit.R
echo "end $(date -u)"
