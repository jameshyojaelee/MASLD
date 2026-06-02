#!/bin/bash
#SBATCH --job-name=net_262_embed
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "262 graph_embedding | Start: $(date)"
python 262_graph_embedding.py
echo "262 graph_embedding | Exit: $? | End: $(date)"
