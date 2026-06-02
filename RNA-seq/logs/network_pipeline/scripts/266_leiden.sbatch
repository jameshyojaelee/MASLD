#!/bin/bash
#SBATCH --job-name=net_266_leiden
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "266 | Start: $(date)"
python 266_*.py
echo "266 | Exit: $? | End: $(date)"
