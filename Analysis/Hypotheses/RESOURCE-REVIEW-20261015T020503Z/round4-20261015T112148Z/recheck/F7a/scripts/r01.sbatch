#!/bin/bash
#SBATCH --job-name=9610
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F7a/logs/9610-%j.out
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F7a
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python $D/scripts/r01_mpra_boot.py $D/out/mpra_$(date -u +%Y%m%dT%H%M%SZ)
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/pip freeze > $D/logs/pip_freeze.txt
