#!/bin/bash
#SBATCH --job-name=netv2_294_portal
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "294 | Start: $(date)"
python 294_export_portal_v2.py
echo "294 | Exit: $? | End: $(date)"
