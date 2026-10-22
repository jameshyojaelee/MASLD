#!/bin/bash
#SBATCH --job-name=H36_build
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H36-20261009T180345Z/logs/slurm-build-%j.out

# H36 step: build. One job, no array. Implements PRESPEC_H36.md.
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=$REPO/Analysis/Hypotheses/H36-20261009T180345Z
SC=$REPO/scripts/analysis/hypotheses_2026q4/H36_splicing
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MASLD_PROJECT_ROOT=$REPO
echo "host $(hostname)  start $(date -u +%FT%TZ)  job ${SLURM_JOB_ID:-NA}"
test -f "$OUT/PRESPEC_H36.md"
test ! -e "$OUT/build"   # never overwrite an earlier run
mkdir -p "$OUT/build"
sha256sum "$OUT/PRESPEC_H36.md" "$SC"/h36_01_build_clusters.R | tee "$OUT/logs/script_sha256_build_${SLURM_JOB_ID:-NA}.txt"
"$RS" "$SC/h36_01_build_clusters.R" "$REPO" "$OUT/build"
echo "end $(date -u +%FT%TZ)"
