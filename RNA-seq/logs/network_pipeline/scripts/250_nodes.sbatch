#!/bin/bash
#SBATCH --job-name=net_250_nodes
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "250 define_network_nodes | Start: $(date)"
Rscript 250_define_network_nodes.R
echo "250 define_network_nodes | Exit: $? | End: $(date)"
