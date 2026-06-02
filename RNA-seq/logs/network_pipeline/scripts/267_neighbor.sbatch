#!/bin/bash
#SBATCH --job-name=net_267_neighbor
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "267 | Start: $(date)"
python 267_*.py
echo "267 | Exit: $? | End: $(date)"
