#!/bin/bash
#SBATCH --job-name=ncRNA_M2_cerna
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=6:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "Module 2: ceRNA Network | Start: $(date)"
Rscript 54_cerna_network.R
echo "Module 2: Exit code: $? | End: $(date)"
