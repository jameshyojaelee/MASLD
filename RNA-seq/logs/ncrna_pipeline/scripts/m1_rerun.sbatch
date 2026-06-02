#!/bin/bash
#SBATCH --job-name=ncRNA_M1_rerun
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "Module 1 RERUN | Start: $(date)"
Rscript 53_ncrna_landscape.R
echo "Module 1 RERUN: Exit code: $? | End: $(date)"
