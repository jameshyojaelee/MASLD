#!/bin/bash
#SBATCH --job-name=N6_statin_rnaseq
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N6-20261009T135702Z/logs/rnaseq_apply-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N6-20261009T135702Z/logs/rnaseq_apply-%j.err
# N6 application step: rank score of the array-derived K=50 gene set in the five RNA-seq cohorts (774 roster units),
# association with histological group and stage, and the refit with the score as a covariate. Single job, no array.
# Usage: SIGNATURE=<signature_K50_all_arrays.tsv> sbatch run_rnaseq_apply.sbatch
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
N6=$ROOT/Analysis/Hypotheses/N6-20261009T135702Z
S=$ROOT/scripts/analysis/hypotheses_2026q4/N6_statin_classifier
ARM=$ROOT/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls
export PATH=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${SIGNATURE:?set SIGNATURE to the signature_K50_all_arrays.tsv path}"
TS=$(date -u +%Y%m%dT%H%M%SZ)
cd "$ROOT"
Rscript "$S/04_rnaseq_apply.R" "$SIGNATURE" "$ARM" "$N6/results/rnaseq_apply-$TS"
echo "done $TS"
