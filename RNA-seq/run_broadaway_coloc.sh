#!/bin/bash
#SBATCH --job-name=liver_broadaway_coloc
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/broadaway_coloc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/broadaway_coloc_%j.err

# Script 35: Broadaway eQTL (N=1,183) COLOC with Ghodsian MASLD GWAS
# Requires bigmem: per-chr eQTL files are ~500MB-1.2GB each
# Estimated runtime: 4-8h with vectorized liftover (was >24h with for-loop)

set -euo pipefail

echo "=== Broadaway COLOC ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

# Setup environment
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/logs

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

echo "Running Script 35..."
Rscript RNA-seq/35_broadaway_coloc.R

echo ""
echo "=== Complete ==="
echo "End: $(date)"
