#!/bin/bash
#SBATCH --job-name=H35_build
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/H35-20261009T180343Z/logs/build-%j.out
#SBATCH --error=Analysis/Hypotheses/H35-20261009T180343Z/logs/build-%j.err
# H35 step 1: person table from the reads extraction and STAR ReadsPerGene.
set -euo pipefail
OUT="${1:?output dir required (must not exist)}"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$PY" -I scripts/analysis/hypotheses_2026q4/H35_editing/01_build_person_table.py "$OUT"
"$PY" -I -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
sha256sum scripts/analysis/hypotheses_2026q4/H35_editing/01_build_person_table.py > "$OUT/script_sha256.txt"
