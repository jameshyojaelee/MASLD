#!/bin/bash
#SBATCH --job-name=H21_sets
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
# H21 step 3: MAGMA competitive gene-set regression for 117 programs and the
# matched random sets, then the difference test and the recomputation check.
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err 03_geneset.sbatch <OUT>
set -euo pipefail
OUT=$1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$ROOT/scripts/analysis/hypotheses_2026q4/H21_program_magma
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
MAGMA=$ROOT/tools/magma/magma
export PYTHONNOUSERSITE=1
cd "$OUT"
mkdir -p geneset geneset_mhc_kept

# Common gene universe: gene IDs only are read from the gene analysis output.
$PY - "$OUT" <<'PYEOF'
import sys, pandas as pd
out = sys.argv[1]
ids = None
for t in ["alt", "pdff", "nafld", "cirrhosis"]:
    g = set(pd.read_csv(f"{out}/gene_analysis/{t}.genes.out", sep=r"\s+", usecols=["GENE"], dtype=str)["GENE"])
    print(t, len(g)); ids = g if ids is None else ids & g
cov = set(pd.read_csv(f"{out}/inputs/gene_covar.tsv", sep="\t", usecols=["GENE"], dtype=str)["GENE"])
mhc = set(pd.read_csv(f"{out}/inputs/mhc_genes.txt", header=None, dtype=str)[0])
with_mhc = sorted(ids & cov); no_mhc = sorted((ids & cov) - mhc)
open(f"{out}/inputs/gene_include.txt", "w").write("\n".join(no_mhc) + "\n")
open(f"{out}/inputs/gene_include_mhc_kept.txt", "w").write("\n".join(with_mhc) + "\n")
print("common genes with covariates:", len(with_mhc), "after MHC exclusion:", len(no_mhc))
PYEOF

for t in alt pdff nafld cirrhosis cirrhosis_mvp bmi; do
  $MAGMA --gene-results gene_analysis/$t.genes.raw \
    --set-annot inputs/all_sets.txt \
    --gene-covar inputs/gene_covar.tsv \
    --model condition-hide=hep_spec,mean_expr \
    --settings gene-include=inputs/gene_include.txt gene-info abbreviate=0 \
    --seed 20261007 \
    --out geneset/$t > logs/geneset_$t.stdout 2>&1
done
for t in alt pdff nafld cirrhosis; do
  $MAGMA --gene-results gene_analysis/$t.genes.raw \
    --set-annot inputs/program_sets.txt \
    --gene-covar inputs/gene_covar.tsv \
    --model condition-hide=hep_spec,mean_expr \
    --settings gene-include=inputs/gene_include_mhc_kept.txt abbreviate=0 \
    --seed 20261007 \
    --out geneset_mhc_kept/$t > logs/geneset_mhc_kept_$t.stdout 2>&1
done

$PY "$HERE/04_difference_and_summary.py" --out "$OUT" > "$OUT/logs/04_summary.stdout"
$PY "$HERE/05_recompute_check.py" --out "$OUT" > "$OUT/results/RECOMPUTE_CHECK.txt"
echo "H21 gene-set analysis done"
