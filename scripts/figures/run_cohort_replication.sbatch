#!/bin/bash
#SBATCH --job-name=cohort_rep
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/output/cohort_replication/slurm_%j.log

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript scripts/figures/fig_cohort_replication.R
