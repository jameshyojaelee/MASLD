#!/bin/bash
#SBATCH --job-name=net_273_idnorm
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "273 ID normalization | Start: $(date)"
python 273_normalize_ids.py
echo "273 ID normalization | Exit: $? | End: $(date)"
