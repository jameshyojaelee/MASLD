#!/bin/bash
#SBATCH --job-name=net_266_annotate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "266 graph_annotation | Start: $(date)"
python 266_graph_annotation.py
echo "266 graph_annotation | Exit: $? | End: $(date)"
