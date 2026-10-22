#!/bin/bash
#SBATCH --job-name=H21_genes
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
# H21 step 2: build inputs, annotate SNPs to genes (35 kb up, 10 kb down),
# run MAGMA SNP-wise mean gene analysis for six GWAS. Single job, no array.
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err 02_gene_analysis.sbatch <OUT>
set -euo pipefail
OUT=$1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$ROOT/scripts/analysis/hypotheses_2026q4/H21_program_magma
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
MAGMA=$ROOT/tools/magma/magma
REF=$ROOT/tools/magma/ref
export PYTHONNOUSERSITE=1
cd "$OUT"

$PY -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
$MAGMA --version > "$OUT/magma_version.txt" 2>&1 || true

$PY "$HERE/01_prepare_inputs.py" --out "$OUT"

$MAGMA --annotate window=35,10 \
  --snp-loc $REF/g1000_eur.bim \
  --gene-loc $REF/NCBI37.3.gene.loc \
  --out "$OUT/gene_analysis/annot_35k10k"

run_one() {
  t=$1
  $MAGMA --bfile $REF/g1000_eur \
    --pval "$OUT/inputs/gwas/$t.pval.tsv" use=SNP,P ncol=N \
    --gene-annot "$OUT/gene_analysis/annot_35k10k.genes.annot" \
    --gene-model snp-wise=mean \
    --seed 20261007 \
    --out "$OUT/gene_analysis/$t" > "$OUT/logs/gene_$t.stdout" 2>&1
}
export -f run_one
export OUT MAGMA REF
printf "%s\n" pdff cirrhosis nafld alt cirrhosis_mvp bmi | xargs -P 3 -I{} bash -c 'run_one {}'
ls -la "$OUT/gene_analysis"
echo "H21 gene analysis done"
