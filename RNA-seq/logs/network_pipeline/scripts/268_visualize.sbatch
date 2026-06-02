#!/bin/bash
#SBATCH --job-name=net_268_visualize
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "268 network_visualization | Start: $(date)"
Rscript 268_network_visualization.R
echo "268 network_visualization | Exit: $? | End: $(date)"
