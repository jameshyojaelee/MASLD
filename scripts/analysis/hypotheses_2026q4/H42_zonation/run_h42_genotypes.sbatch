#!/bin/bash
#SBATCH --job-name=H42_geno
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/H42-20261009T180511Z/logs/geno-%j.out
#SBATCH --error=Analysis/Hypotheses/H42-20261009T180511Z/logs/geno-%j.err
# H42 step 1: person-level PNPLA3 rs738409 genotype calls from the finished reads extraction.
set -euo pipefail
H=Analysis/Hypotheses/H42-20261009T180511Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
OUT=$H/results/geno-$(date -u +%Y%m%dT%H%M%SZ)
"$PY" -I scripts/analysis/hypotheses_2026q4/H42_zonation/01_call_genotypes.py \
  Analysis/Hypotheses/READS-20261009T141017Z \
  figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv \
  "$OUT"
"$PY" -m pip freeze > "$OUT/pip_freeze.txt"
echo "OUT=$OUT"
