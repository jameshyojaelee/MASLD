#!/bin/bash
#SBATCH --job-name=netv2_284_dlr
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "284 | Start: $(date)"
Rscript 284_lr_differential.R
echo "284 | Exit: $? | End: $(date)"
