#!/bin/bash
#SBATCH --job-name=8035
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -o /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H63/8035_%j.out
W=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H63
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 TMPDIR=$W/tmp
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python $W/count_detect.py $W
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m pip freeze > $W/pip_freeze.txt
