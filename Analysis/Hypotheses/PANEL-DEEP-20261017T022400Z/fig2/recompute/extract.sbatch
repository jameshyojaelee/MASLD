#!/bin/bash
#SBATCH --job-name=9820
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -p cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig2/recompute/extract_%j.out
export PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig2/recompute/extract_upstream.py
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/pip freeze > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig2/recompute/pip_freeze.txt
