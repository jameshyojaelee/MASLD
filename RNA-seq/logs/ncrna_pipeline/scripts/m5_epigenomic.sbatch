#!/bin/bash
#SBATCH --job-name=ncRNA_M5_epigenomic
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "Module 5: Epigenomic Regulation | Start: $(date)"
Rscript 56_ncrna_epigenomic.R
echo "Module 5: Exit code: $? | End: $(date)"
