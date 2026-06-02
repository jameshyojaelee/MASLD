#!/bin/bash
#SBATCH --job-name=netv2_281_f2_delta
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "281 | Start: $(date)"
Rscript 281_f2_delta_edges.R
echo "281 | Exit: $? | End: $(date)"
