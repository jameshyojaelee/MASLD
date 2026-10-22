#!/bin/bash
#SBATCH --job-name=9812
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1/slurm-9812-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1/slurm-9812-%j.err
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
PYTHONNOUSERSITE=1 $PY Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1/recompute_1c_untested.py
