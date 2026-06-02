#!/bin/bash
#SBATCH --job-name=net_261_threshold
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "261 edge_thresholding | Start: $(date)"
python 261_edge_thresholding.py
echo "261 edge_thresholding | Exit: $? | End: $(date)"
