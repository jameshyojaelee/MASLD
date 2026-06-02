#!/bin/bash
#SBATCH --job-name=pathway_programs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pathway_programs_%j.log

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
source ~/.bashrc
micromamba activate rnaseq

Rscript scripts/figures/figS_pathway_programs.R
