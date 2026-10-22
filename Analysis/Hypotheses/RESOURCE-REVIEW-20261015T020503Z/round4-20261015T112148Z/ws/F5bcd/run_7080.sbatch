#!/bin/bash
#SBATCH --job-name=7080
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 2
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5bcd/logs/7080_%j.out
WS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5bcd
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 TMPDIR=$WS/tmp
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY $WS/fig5d_locus_background.py $WS/fig5d_locus_background_rows.tsv
$PY -m pip freeze > $WS/pip_freeze_7080.txt
