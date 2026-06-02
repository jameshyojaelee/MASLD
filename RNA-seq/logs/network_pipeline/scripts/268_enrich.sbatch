#!/bin/bash
#SBATCH --job-name=net_268_enrich
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "268 | Start: $(date)"
Rscript 268_*.R
echo "268 | Exit: $? | End: $(date)"
