#!/bin/bash
#SBATCH --job-name=stg171_contrastive
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/171_contrastive_%j.out
#SBATCH --error=logs/171_contrastive_%j.err

set -euo pipefail

echo "=== Script 171: Contrastive Multi-Source Gene Embeddings ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Date:   $(date)"
echo ""

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

# Create output dir
mkdir -p ../results/novel_ml/contrastive_embeddings

# Create logs dir if needed
mkdir -p logs

python 171_contrastive_gene_embeddings.py

echo ""
echo "=== Done: $(date) ==="
