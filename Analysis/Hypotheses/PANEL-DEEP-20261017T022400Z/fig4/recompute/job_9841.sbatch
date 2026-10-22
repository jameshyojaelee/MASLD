#!/bin/bash
#SBATCH --job-name=9841
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute/job_9841_%j.log
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONNOUSERSITE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python r4gh_recompute.py
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m pip freeze > pip_freeze_9841.txt
