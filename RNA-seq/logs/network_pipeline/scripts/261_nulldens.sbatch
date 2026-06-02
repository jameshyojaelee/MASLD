#!/bin/bash
#SBATCH --job-name=net_261_nulldens
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "261 | Start: $(date)"
python 261_null_density_estimation.py
echo "261 | Exit: $? | End: $(date)"
