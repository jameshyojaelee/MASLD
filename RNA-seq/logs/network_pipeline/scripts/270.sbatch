#!/bin/bash
#SBATCH --job-name=net_270_atlas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "270 270_atlas_integration.R | Start: $(date)"
Rscript 270_atlas_integration.R
echo "270 270_atlas_integration.R | Exit: $? | End: $(date)"
