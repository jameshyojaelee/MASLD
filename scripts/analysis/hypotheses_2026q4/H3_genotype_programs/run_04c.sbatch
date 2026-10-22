#!/bin/bash
#SBATCH --job-name=H3_sig_comb
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/H3-20261009T180331Z/logs/%x-%j.out
set -euo pipefail
ROOT=Analysis/Hypotheses/H3-20261009T180331Z
PN=$(cat $ROOT/results/LATEST_SIGNATURE_pnpla3.txt)
HS=$(cat $ROOT/results/LATEST_SIGNATURE_hsd.txt)
OUT=$ROOT/results/signature-combined-$(date -u +%Y%m%dT%H%M%SZ)
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -I scripts/analysis/hypotheses_2026q4/H3_genotype_programs/04c_combine_signature.py "$PN" "$HS" "$OUT"
sha256sum scripts/analysis/hypotheses_2026q4/H3_genotype_programs/04c_combine_signature.py > "$OUT/script_sha256.txt"
echo "$OUT" > $ROOT/results/LATEST_SIGNATURE.txt
