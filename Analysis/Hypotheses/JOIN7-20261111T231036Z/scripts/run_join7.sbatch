#!/bin/bash
#SBATCH --job-name=join7
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/JOIN7-20261111T231036Z/logs/join7_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/JOIN7-20261111T231036Z/logs/join7_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$ROOT/Analysis/Hypotheses/JOIN7-20261111T231036Z
OUT=$HERE/out_${SLURM_JOB_ID}
test ! -e $OUT || { echo "output exists: $OUT"; exit 1; }
mkdir -p $OUT
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=20261111
cd $ROOT
{ $PY --version; hostname; date -u; } > $OUT/versions.txt
$PY -m pip freeze > $OUT/pip_freeze.txt
sha256sum GWAS/MR_Data/GCST90267352_PDFF_Pazoki2022.tsv.gz GWAS/MR_Data/Ghouse_Cirrhosis/GCST90319877_harmonised_hg38.tsv.gz \
  GWAS/MR_Data/PanUKBB/biomarkers-30620-ALT-both_sexes-irnt.tsv.bgz GWAS/MR_Data/MVP/raw/labs/MVP_R4.1000G_AGR.ALT_Mean_INT.EUR.GIA.dbGaP.txt.gz \
  tools/pops/example/data/utils/gene_annot_jun10.txt Analysis/Hypotheses/GSET8-20261110T235125Z/out_23188024/pdff.genes.out \
  Analysis/Hypotheses/GSET8-20261110T235125Z/out_23188024/cirrhosis.genes.out \
  Analysis/Hypotheses/DIR6-20261111T024808Z/out_23210381/dir6_gene_classes.tsv Analysis/Hypotheses/DIR6-20261111T024808Z/out_23210381/calibration.tsv \
  Analysis/Hypotheses/TRIAL1-20261110T205624Z/trial_outcomes_v2.tsv Analysis/Hypotheses/TRIAL1-20261110T205624Z/excluded_v2.tsv \
  Analysis/Hypotheses/TRIAL1-20261110T205624Z/CHANGES_v2.tsv figures/main/fig6_gene_catalog/panels/fig6_drug_target_calibration_values.tsv \
  RNA-seq/results/evidence_passports/candidates/gene-catalog-r8-anchor-disease-sidecar-20261025T131731Z/catalog_r8_anchor_disease_liability.tsv \
  > $OUT/input_sha256.txt
{ find $HERE/scripts -maxdepth 1 -type f | sort | xargs sha256sum; sha256sum $HERE/PRESPEC_JOIN7.md; } > $OUT/script_sha256.txt
$PY $HERE/scripts/10_extract.py $ROOT $OUT
$PY $HERE/scripts/20_join.py $ROOT $OUT
date -u >> $OUT/versions.txt
