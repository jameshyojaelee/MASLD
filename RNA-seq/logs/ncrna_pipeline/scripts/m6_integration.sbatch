#!/bin/bash
#SBATCH --job-name=ncRNA_M6_integration
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "Module 6: Atlas Integration | Start: $(date)"
Rscript 57_ncrna_atlas_integration.R
echo "Module 6: Exit code: $? | End: $(date)"
