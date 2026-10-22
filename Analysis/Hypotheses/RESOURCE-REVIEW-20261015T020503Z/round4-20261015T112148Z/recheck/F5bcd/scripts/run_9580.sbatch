#!/bin/bash
#SBATCH --job-name=9580
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F5bcd/logs/9580_%j.out
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F5bcd/scripts/r11_02_background_indep.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F5bcd/out
