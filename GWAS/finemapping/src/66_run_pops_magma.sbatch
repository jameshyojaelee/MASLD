#!/bin/bash
#SBATCH --job-name=magma
#SBATCH --partition=cpu
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/66_magma_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/66_magma_%j.err

# PoPS Axis-2 step 1: ENSGID-keyed MAGMA gene analysis on the primary MASLD GWAS
# (MVP_NAFLD_EUR). ADDITIVE. Reuses tools/magma binary + g1000_eur bfile/snploc.
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p GWAS/finemapping/results/seqfunc/logs

MAGMA_BIN=tools/magma/magma
REF=tools/magma/ref
WORK=tools/magma/work
OUT=$WORK/pops
mkdir -p $OUT

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

echo "[magma] step 1: prep ENSGID gene.loc + MVP_NAFLD_EUR MAGMA input"
python GWAS/finemapping/src/66_pops_magma_prep.py

echo "[magma] step 2: MAGMA --annotate (SNP->gene, ENSGID, 0kb window)"
$MAGMA_BIN \
  --annotate \
  --snp-loc $WORK/g1000_eur.snploc \
  --gene-loc $OUT/gene_annot_jun10.ensgid.gene.loc \
  --out $OUT/magma_annot_ensgid

echo "[magma] step 3: MAGMA gene analysis on MVP_NAFLD_EUR"
$MAGMA_BIN \
  --bfile $REF/g1000_eur \
  --pval $OUT/MVP_NAFLD_EUR.magma_input.tsv use=SNP,P ncol=N \
  --gene-annot $OUT/magma_annot_ensgid.genes.annot \
  --out $OUT/MVP_NAFLD_EUR

echo "[magma] done."
echo "  genes.out: $OUT/MVP_NAFLD_EUR.genes.out"
head -3 $OUT/MVP_NAFLD_EUR.genes.out || true
wc -l $OUT/MVP_NAFLD_EUR.genes.out || true
