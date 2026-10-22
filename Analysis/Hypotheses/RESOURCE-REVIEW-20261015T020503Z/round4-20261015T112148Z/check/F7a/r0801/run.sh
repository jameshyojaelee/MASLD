#!/bin/bash
#SBATCH -J 9110
#SBATCH -p cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -o /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F7a/r0801/slurm_%j.out
export PYTHONDONTWRITEBYTECODE=1 TMPDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F7a/tmp
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F7a/r0801/mpra7.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F7a/r0801/mpra7_result.tsv
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/pip freeze > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F7a/r0801/pip_freeze.txt
