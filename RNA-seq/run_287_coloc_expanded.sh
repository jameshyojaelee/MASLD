#!/bin/bash
#SBATCH --job-name=net_287_dcoloc_expanded
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=logs/net_287_dcoloc_expanded_%j.out
#SBATCH --error=logs/net_287_dcoloc_expanded_%j.err

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

echo "287 expanded D-COLOC edges | Start: $(date)"
Rscript 287_coloc_locus_edges_expanded.R
echo "287 expanded D-COLOC edges | Exit: $? | End: $(date)"
