#!/bin/bash
#SBATCH --job-name=broadaway_multi_coloc
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/broadaway_%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/broadaway_%x_%j.err

# Script 35g: Broadaway eQTL (N=1,183) COLOC with parameterized UKBB GWAS
# Requires bigmem: per-chr eQTL files are ~500MB-1.2GB each
#
# Usage:
#   GWAS_NAME=UKBB_AST sbatch RNA-seq/run_broadaway_coloc_multi_gwas.sh
#   GWAS_NAME=UKBB_GGT sbatch RNA-seq/run_broadaway_coloc_multi_gwas.sh
#   GWAS_NAME=UKBB_ALT sbatch RNA-seq/run_broadaway_coloc_multi_gwas.sh  (equivalent to 35b)

set -euo pipefail

GWAS_NAME="${GWAS_NAME:-UKBB_AST}"
export GWAS_NAME

echo "=== Broadaway x ${GWAS_NAME} COLOC ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

# Setup environment
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs "RNA-seq/results/causal_inference/broadaway_$(echo "$GWAS_NAME" | tr '[:upper:]' '[:lower:]')"

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

echo "Running Script 35g with GWAS_NAME=${GWAS_NAME}..."
Rscript RNA-seq/35g_broadaway_coloc_multi_gwas.R

echo ""
echo "=== Complete ==="
echo "End: $(date)"
