#!/bin/bash
#SBATCH --job-name=netv2_285_dcerna
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "285 | Start: $(date)"
Rscript 285_cerna_masld_expressed.R
echo "285 | Exit: $? | End: $(date)"
