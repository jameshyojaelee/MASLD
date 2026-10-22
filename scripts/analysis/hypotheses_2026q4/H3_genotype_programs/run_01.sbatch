#!/bin/bash
#SBATCH --job-name=H3_geno
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/H3-20261009T180331Z/logs/%x-%j.out
set -euo pipefail
ROOT=Analysis/Hypotheses/H3-20261009T180331Z
OUT=$ROOT/results/geno-$(date -u +%Y%m%dT%H%M%SZ)
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY -I scripts/analysis/hypotheses_2026q4/H3_genotype_programs/01_person_genotypes.py "$OUT"
$PY -m pip freeze > "$OUT/pip_freeze.txt"
sha256sum scripts/analysis/hypotheses_2026q4/H3_genotype_programs/01_person_genotypes.py > "$OUT/script_sha256.txt"
echo "$OUT" > $ROOT/results/LATEST_GENO.txt
