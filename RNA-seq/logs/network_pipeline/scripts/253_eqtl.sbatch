#!/bin/bash
#SBATCH --job-name=net_253_eqtl
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "253 eqtl_edges | Start: $(date)"
Rscript 253_eqtl_edges.R
echo "253 eqtl_edges | Exit: $? | End: $(date)"
