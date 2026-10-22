#!/bin/bash
#SBATCH --job-name=H3_sig_hsd
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/H3-20261009T180331Z/logs/%x-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1
ROOT=Analysis/Hypotheses/H3-20261009T180331Z
CALLS=$(cat $ROOT/results/LATEST_CALLS.txt)
SCORES=$(cat $ROOT/results/LATEST_SCORES.txt)
OUT=$ROOT/results/signature-hsd-$(date -u +%Y%m%dT%H%M%SZ)
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
$RS scripts/analysis/hypotheses_2026q4/H3_genotype_programs/04b_carrier_signature_model.R "$CALLS" "$SCORES" "$OUT" hsd
sha256sum scripts/analysis/hypotheses_2026q4/H3_genotype_programs/04b_carrier_signature_model.R > "$OUT/script_sha256.txt"
echo "$OUT" > $ROOT/results/LATEST_SIGNATURE_hsd.txt
