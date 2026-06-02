#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --job-name=fig4_test
#SBATCH --output=logs/fig4_test_%j.out
#SBATCH --error=logs/fig4_test_%j.err

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript scripts/figures/fig5_causal_architecture.R
