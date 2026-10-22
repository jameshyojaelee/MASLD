#!/bin/bash
#SBATCH --job-name=9073
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72/9073_%j.out
W=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72
export TMPDIR=$W/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python $W/mde.py
sha256sum /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/MR_Data/PanUKBB/sex_stratified/30*_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz $W/mde.py $W/mde.out $W/best_enzyme_row.tsv > $W/mde_sha256.txt
