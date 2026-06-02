#!/bin/bash
#SBATCH --job-name=netv2_281b_f2
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "281b | Start: $(date)"
Rscript 281b_f2_delta_edges_pragmatic.R
echo "281b | Exit: $? | End: $(date)"
