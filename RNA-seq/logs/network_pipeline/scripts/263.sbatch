#!/bin/bash
#SBATCH --job-name=net_263_noisyor
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00


eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "263 263_noisy_or_composite.py | Start: $(date)"
python 263_noisy_or_composite.py
echo "263 263_noisy_or_composite.py | Exit: $? | End: $(date)"
