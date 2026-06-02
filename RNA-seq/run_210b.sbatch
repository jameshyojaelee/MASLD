#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=210b_prog_plots
#SBATCH --output=logs/210b_%j.out
#SBATCH --error=logs/210b_%j.err

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
Rscript 210b_progression_coloc_plots.R
