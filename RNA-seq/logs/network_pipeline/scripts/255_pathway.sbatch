#!/bin/bash
#SBATCH --job-name=net_255_pathway
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "255 pathway_edges | Start: $(date)"
Rscript 255_pathway_edges.R
echo "255 pathway_edges | Exit: $? | End: $(date)"
