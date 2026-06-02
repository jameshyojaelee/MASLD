#!/bin/bash
#SBATCH --job-name=net_267_prioritize
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "267 node_prioritization | Start: $(date)"
python 267_node_prioritization.py
echo "267 node_prioritization | Exit: $? | End: $(date)"
