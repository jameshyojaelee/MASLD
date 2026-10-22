#!/bin/bash
#SBATCH --job-name=9834
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -p cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/logs/9834_%j.out
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/code/19_lncrna_matched.R
