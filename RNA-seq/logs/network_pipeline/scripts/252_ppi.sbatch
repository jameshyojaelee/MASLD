#!/bin/bash
#SBATCH --job-name=net_252_ppi
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "252 ppi_edges | Start: $(date)"
Rscript 252_ppi_edges.R
echo "252 ppi_edges | Exit: $? | End: $(date)"
