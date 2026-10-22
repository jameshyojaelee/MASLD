#!/bin/bash
#SBATCH --job-name=9540
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F3cdg/run-20261015T173135Z/slurm-%j.out
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
E=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
O=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F3cdg/run-20261015T173135Z
$E/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F3cdg/recheck_refits.R $O
PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 $E/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F3cdg/recheck_tables.py $O
PYTHONNOUSERSITE=1 $E/python -m pip freeze > $O/pip_freeze.txt
