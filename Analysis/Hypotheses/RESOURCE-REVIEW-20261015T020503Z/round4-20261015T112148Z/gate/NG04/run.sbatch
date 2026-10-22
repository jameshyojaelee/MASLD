#!/bin/bash
#SBATCH -J 8050
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH --time 00:30:00
#SBATCH -p cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -o /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG04/8050_%j.log
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG04
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 TMPDIR=$D/tmp MPLCONFIGDIR=$D/tmp
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY $D/gate_count.py
$PY -m pip freeze > $D/pip_freeze.txt
