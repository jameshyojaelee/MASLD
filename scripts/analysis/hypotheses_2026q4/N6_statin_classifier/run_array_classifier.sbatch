#!/bin/bash
#SBATCH --job-name=N6_statin_array
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N6-20261009T135702Z/logs/array_classifier-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N6-20261009T135702Z/logs/array_classifier-%j.err
# N6 steps 2 and 3: gene-level matrix from GSE130991, then repeated CV, permutation nulls,
# confound checks, positive control, gate G2. Single job, no array. Output dirs are new and timestamped.
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
N6=$ROOT/Analysis/Hypotheses/N6-20261009T135702Z
S=$ROOT/scripts/analysis/hypotheses_2026q4/N6_statin_classifier
export PATH=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TS=$(date -u +%Y%m%dT%H%M%SZ)
cd "$ROOT"
Rscript "$S/02_array_prepare.R" "$N6/downloads_untrusted/GSE130991_series_matrix.txt.gz" \
  "$N6/downloads_untrusted/GPL20265_platform_table.txt" "${SAMPLES:-$N6/results/samples.tsv}" "$N6/results/prep-$TS"
Rscript "$S/03_array_classifier.R" "$N6/results/prep-$TS/gene_matrix.rds" "${SAMPLES:-$N6/results/samples.tsv}" \
  "$N6/results/classifier-$TS"
echo "done $TS"
