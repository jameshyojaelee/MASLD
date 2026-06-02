#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --job-name=fig4
#SBATCH --output=scripts/figures/logs/fig4_%j.out
#SBATCH --error=scripts/figures/logs/fig4_%j.err

source ~/.bashrc
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript scripts/figures/fig4_causal_architecture.R
