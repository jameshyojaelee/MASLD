#!/bin/bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH -o /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F4abcd/logs/%x_%j.out
#SBATCH -e /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F4abcd/logs/%x_%j.err
CK=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F4abcd
export TMPDIR=$CK/tmp PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY -m pip freeze > $CK/out/content/pip_freeze.txt 2>/dev/null
/usr/bin/time -v $PY $CK/scripts/check_content_bootstrap.py $CK/out/content
