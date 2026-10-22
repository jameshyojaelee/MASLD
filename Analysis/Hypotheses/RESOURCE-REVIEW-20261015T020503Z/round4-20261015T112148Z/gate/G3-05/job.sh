#!/bin/bash
#SBATCH -J 8230
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH --time 00:30:00
#SBATCH -p cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -o /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/G3-05/job_%j.out
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/G3-05
export PYTHONDONTWRITEBYTECODE=1 TMPDIR=$PWD
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python count_coverage.py
