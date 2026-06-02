#!/bin/bash
#SBATCH --job-name=net_270_atlas
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "270 | Start: $(date)"
Rscript 270_*.R
echo "270 | Exit: $? | End: $(date)"
