#!/bin/bash
#SBATCH --job-name=net_264_diagnost
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "264 | Start: $(date)"
python 264_*.py
echo "264 | Exit: $? | End: $(date)"
