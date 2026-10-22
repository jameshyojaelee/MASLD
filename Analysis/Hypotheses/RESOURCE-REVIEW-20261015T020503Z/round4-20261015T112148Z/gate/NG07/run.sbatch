#!/bin/bash
#SBATCH --job-name=8080
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG07/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG07/logs/%x_%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 TMPDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG07/tmp
RP=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
LE=/gpfs/commons/home/jameslee/micromamba/envs/ldsc/bin
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG07/ldsc
$RP /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG07/02_panukbb_to_ldsc_input.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/MR_Data/PanUKBB/biomarkers-30730-GGT-both_sexes-irnt.tsv.bgz /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/ldsc_ref/hm3_chrbp_to_rsid.tsv /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/ldsc_ref/hm3_merge_alleles.txt 400751 /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG07/ldsc/ggt_src.tsv.gz
$LE/python $LE/munge_sumstats.py --sumstats /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG07/ldsc/ggt_src.tsv.gz --out /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG07/ldsc/ggt --merge-alleles /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/ldsc_ref/hm3_merge_alleles.txt --N-col N --a1 A1 --a2 A2 --snp SNP --p P --signed-sumstats BETA,0 --chunksize 500000
$RP /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG07/gate.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG07 /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H22-20261009T015913Z /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H22-20261009T015913Z/ldsc/sumstats/alt.sumstats.gz /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG07/ldsc/ggt.sumstats.gz
$RP -m pip freeze > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG07/pip_freeze_rnaseq.txt
