#!/bin/bash
#SBATCH --job-name=9842
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute/job_9842_%j.log
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript r4ce_refit.R
