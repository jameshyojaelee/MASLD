#!/bin/bash
#SBATCH --job-name=net_286_dstable
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=logs/net_286_dstable_%j.out
#SBATCH --error=logs/net_286_dstable_%j.err

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

echo "286 stable co-expression edges | Start: $(date)"
Rscript 286_stable_coexpr_edges.R
echo "286 stable co-expression edges | Exit: $? | End: $(date)"
