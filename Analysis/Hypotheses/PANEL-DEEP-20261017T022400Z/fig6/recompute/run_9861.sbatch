#!/bin/bash
#SBATCH --job-name=9861
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig6/recompute/sensitivity_9861_%j.out
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig6/recompute/sensitivity_fig6.py
