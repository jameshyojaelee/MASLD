#!/bin/bash
#SBATCH --job-name=netv2_292_f2_leiden
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "292 | Start: $(date)"
python 292_f2_weighted_leiden.py
echo "292 | Exit: $? | End: $(date)"
