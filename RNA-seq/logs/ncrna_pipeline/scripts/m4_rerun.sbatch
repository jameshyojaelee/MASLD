#!/bin/bash
#SBATCH --job-name=ncRNA_M4_rerun
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "Module 4 RERUN (W2 Ensembl->MGI fix) | Start: $(date)"
Rscript 55_ncrna_conservation.R
echo "Module 4 RERUN: Exit code: $? | End: $(date)"
