#!/bin/bash
#SBATCH --job-name=net_256_drug
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "256 drug_gene_edges | Start: $(date)"
Rscript 256_drug_gene_edges.R
echo "256 drug_gene_edges | Exit: $? | End: $(date)"
