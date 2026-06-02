#!/bin/bash
#SBATCH --job-name=liver_sceqtl_ukbb
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/sceqtl_coloc_ukbb_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/sceqtl_coloc_ukbb_%j.err

# Script 36b: UKBB ALT COLOC with sc-eQTL (all eGenes, no GWAS gate)
# 4 cell types × ~9,000 eGenes each; estimated 4-6h

set -euo pipefail

echo "=== sc-eQTL COLOC with UKBB ALT GWAS ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs

echo "Running Script 36b..."
Rscript RNA-seq/36b_sceqtl_coloc_ukbb.R

echo ""
echo "=== Complete ==="
echo "End: $(date)"
