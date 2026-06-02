#!/bin/bash
#SBATCH --job-name=net_256_spat_cos
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "256 256_spatial_cosmos_xspecies_edges.R | Start: $(date)"
Rscript 256_spatial_cosmos_xspecies_edges.R
echo "256 256_spatial_cosmos_xspecies_edges.R | Exit: $? | End: $(date)"
