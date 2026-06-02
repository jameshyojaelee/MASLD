#!/bin/bash
#SBATCH --job-name=net_262_postcomp
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "262 | Start: $(date)"
python 262_compute_posteriors.py
echo "262 | Exit: $? | End: $(date)"
