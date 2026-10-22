#!/bin/bash
#SBATCH --job-name=9530
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F2d/run-20261015T174052Z/logs/9530-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F2d/run-20261015T174052Z/logs/9530-%j.err
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F2d/run-20261015T174052Z/scripts/rederive_recheck.py
$PY -m pip freeze > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F2d/run-20261015T174052Z/pip_freeze.txt
