#!/bin/bash
#SBATCH --job-name=9080
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 2
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH -o /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F5bcd/9080_%j.log
W=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F5bcd
export PYTHONDONTWRITEBYTECODE=1 TMPDIR=$W/tmp
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python $W/bg_check.py $W/bg_rows_check.tsv
