#!/bin/bash
#SBATCH --job-name=9072
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72/9072_%j.out
set -uo pipefail
R=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
W=$R/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72
export TMPDIR=$W/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
cd $R
( sha256sum \
  GWAS/MR_Data/PanUKBB/sex_stratified/30620_irnt.gwas.imputed_v3.female.varorder.tsv.bgz \
  GWAS/MR_Data/PanUKBB/sex_stratified/30620_irnt.gwas.imputed_v3.male.varorder.tsv.bgz \
  GWAS/MR_Data/PanUKBB/sex_stratified/30650_irnt.gwas.imputed_v3.female.varorder.tsv.bgz \
  GWAS/MR_Data/PanUKBB/sex_stratified/30650_irnt.gwas.imputed_v3.male.varorder.tsv.bgz \
  GWAS/MR_Data/PanUKBB/sex_stratified/30730_irnt.gwas.imputed_v3.female.varorder.tsv.bgz \
  GWAS/MR_Data/PanUKBB/sex_stratified/30730_irnt.gwas.imputed_v3.male.varorder.tsv.bgz \
  data/public_r4/neale_30880/30880_irnt.gwas.imputed_v3.female.varorder.tsv.bgz \
  data/public_r4/neale_30880/30880_irnt.gwas.imputed_v3.male.varorder.tsv.bgz \
  GWAS/finemapping/data/sumstats/UKBB_ALT_reformatted_hg19.tsv \
  GWAS/finemapping/data/sumstats/UKBB_AST_reformatted_hg19.tsv \
  GWAS/finemapping/data/sumstats/UKBB_GGT_reformatted_hg19.tsv \
  GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv \
  figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet \
  Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L06/ukbb_alt_clumped500kb.tsv \
  Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L06/ukbb_ast_clumped500kb.tsv \
  Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L06/ukbb_ggt_clumped500kb.tsv \
  Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L06-02/clump.py \
  GWAS/finemapping/config/gwas_registry.tsv \
  Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG11/positive_control.txt \
  Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG11/counts.txt \
  Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG11/pc.py \
  > $W/input_sha256.txt; stat -c '%s %n' $(awk '{print $2}' $W/input_sha256.txt) > $W/input_bytes.txt ) &
$PY $W/count_n.py
wait
sha256sum $W/count_n.py $W/count_n.out >> $W/input_sha256.txt
$PY -m pip freeze > $W/pip_freeze.txt 2>/dev/null
