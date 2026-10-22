#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3841
B=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export TMPDIR=$J/tmp; mkdir -p $TMPDIR
cd $J
date
zcat $B/data/external/currin_2025_caqtl/hg19ToHg38.over.chain.gz > hg19ToHg38.over.chain
for s in UKBB_GGT MVP_Platelet_EUR; do
  awk -F'\t' -v o1="${s}_chr1.tsv" -v o19="${s}_chr19.tsv" 'NR>1 && $1=="1"{print > o1} NR>1 && $1=="19"{print > o19}' $B/GWAS/finemapping/data/sumstats/${s}_reformatted_hg19.tsv
  wc -l ${s}_chr1.tsv ${s}_chr19.tsv
done
date
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript pilot_caqtl_coloc.R
date
