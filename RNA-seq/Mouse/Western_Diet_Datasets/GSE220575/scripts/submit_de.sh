#!/bin/bash
#SBATCH --job-name=GSE220575_DE
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets/GSE220575/logs/de_analysis_%j.log

source ~/.bashrc
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets/GSE220575

Rscript scripts/run_de_analysis.R
