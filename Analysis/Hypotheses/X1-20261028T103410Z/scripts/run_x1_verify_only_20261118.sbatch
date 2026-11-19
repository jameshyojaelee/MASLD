#!/bin/bash
#SBATCH --job-name=x1_verify_261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/X1-20261028T103410Z
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
CAT=$ROOT/RNA-seq/results/evidence_passports/candidates/gene-catalog-unified-20261118/catalog_annotated_gene_table_unified_20261118.tsv.gz
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
OUT=out_20261118_unified_cat

export PYTHONNOUSERSITE=1
"$PY" -I scripts/verify_x1_unified_cat_20261118.py "$ROOT/Analysis/Hypotheses/X1-20261028T103410Z/$OUT" "$CAT"
"$PY" -m pip freeze > "$OUT/pip_freeze.txt"
echo "X1_VERIFY_ONLY_DONE"
