#!/bin/bash
#SBATCH --job-name=H3_calls
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/H3-20261009T180331Z/logs/%x-%j.out
set -euo pipefail
ROOT=Analysis/Hypotheses/H3-20261009T180331Z
GENO=$(cat $ROOT/results/LATEST_GENO.txt)
OUT=$ROOT/results/calls-$(date -u +%Y%m%dT%H%M%SZ)
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY -I scripts/analysis/hypotheses_2026q4/H3_genotype_programs/01b_amended_calls.py "$GENO" "$OUT"
$PY -m pip freeze > "$OUT/pip_freeze.txt"
sha256sum scripts/analysis/hypotheses_2026q4/H3_genotype_programs/01b_amended_calls.py > "$OUT/script_sha256.txt"
echo "$OUT" > $ROOT/results/LATEST_CALLS.txt
