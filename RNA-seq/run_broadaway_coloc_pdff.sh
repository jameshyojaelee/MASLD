#!/bin/bash
#SBATCH --job-name=broadaway_pdff_coloc
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/broadaway_pdff_coloc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/broadaway_pdff_coloc_%j.err

# Script 35h: Broadaway eQTL (N=1,183) COLOC with PDFF GWAS (N=33,588)
# Both datasets are GRCh37 — no liftover needed (faster than AST/GGT)
# Requires bigmem: per-chr eQTL files are ~500MB-1.2GB each

set -euo pipefail

echo "=== Broadaway x PDFF COLOC ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

# Setup environment
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs RNA-seq/results/causal_inference/broadaway_pdff

echo "Running Script 35h (no liftover — hg19 native merge)..."
Rscript RNA-seq/35h_broadaway_coloc_pdff.R

echo ""
echo "=== Complete ==="
echo "End: $(date)"
