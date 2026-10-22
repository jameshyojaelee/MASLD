#!/bin/bash
#SBATCH --job-name=H9_locus_classes
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H9-20261009T140001Z/logs/slurm-%j.out

# H9: published MASLD locus classes (Jamialahmadi 2024) versus the 117 Hotspot programs.
# Implements PRESPEC_H9.md. One job, no arrays. Writes to a new results/run-<UTC> directory.
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
H9=$REPO/Analysis/Hypotheses/H9-20261009T140001Z
SC=$REPO/scripts/analysis/hypotheses_2026q4/H9_locus_classes
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "host $(hostname)  start $(date -u +%FT%TZ)  job ${SLURM_JOB_ID:-NA}"
test -f "$H9/PRESPEC_H9.md"
sha256sum "$H9/PRESPEC_H9.md" "$SC/h9_locus_classes.R" "$SC/run_h9.sbatch" \
  | tee "$H9/logs/script_sha256_${SLURM_JOB_ID:-NA}.txt"
"$RS" --vanilla "$SC/h9_locus_classes.R"
echo "end $(date -u +%FT%TZ)"
