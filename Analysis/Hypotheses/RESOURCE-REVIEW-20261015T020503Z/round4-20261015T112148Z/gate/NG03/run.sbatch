#!/bin/bash
#SBATCH -J 8040
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH --time 00:30:00
#SBATCH -p cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -o /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG03/8040-%j.out
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG03
export TMPDIR=$D/tmp PYTHONDONTWRITEBYTECODE=1
P=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$P -m pip freeze > $D/pip_freeze.txt 2>/dev/null
$P $D/gate_ng03.py
