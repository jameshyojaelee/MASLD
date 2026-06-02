#!/bin/bash
#SBATCH --job-name=net_260_permnull
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "260 | Start: $(date)"
python 260_permutation_nulls.py
echo "260 | Exit: $? | End: $(date)"
