#!/bin/bash
#SBATCH --job-name=net_254_genetic
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "254 254_genetic_coloc_edges.R | Start: $(date)"
Rscript 254_genetic_coloc_edges.R
echo "254 254_genetic_coloc_edges.R | Exit: $? | End: $(date)"
