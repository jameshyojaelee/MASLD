#!/bin/bash
#SBATCH --job-name=broadaway_ukbb_coloc
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/broadaway_ukbb_coloc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/broadaway_ukbb_coloc_%j.err

# Script 35b: Broadaway eQTL (N=1,183) COLOC with UKBB ALT GWAS (N=343,850)
# Requires bigmem: per-chr eQTL files are ~500MB-1.2GB each
# Tests ALL 6,564 eGenes (no DEG filter). Estimated runtime: 4-8h.

set -euo pipefail

echo "=== Broadaway x UKBB ALT COLOC ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

# Setup environment
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs RNA-seq/results/causal_inference/broadaway_ukbb

# Download chain file if not present
CHAIN_GZ="data/broadaway_eqtl/hg19ToHg38.over.chain.gz"
CHAIN="data/broadaway_eqtl/hg19ToHg38.over.chain"
if [ ! -f "$CHAIN" ]; then
  if [ ! -f "$CHAIN_GZ" ]; then
    echo "Downloading hg19ToHg38 chain file..."
    wget -q -O "$CHAIN_GZ" \
      "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
  fi
  echo "Decompressing chain file..."
  gunzip -k "$CHAIN_GZ"
fi

echo "Running Script 35b..."
Rscript RNA-seq/35b_broadaway_coloc_ukbb.R

echo ""
echo "=== Complete ==="
echo "End: $(date)"
