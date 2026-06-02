#!/bin/bash
#SBATCH --job-name=script61_features
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=logs/script61_%j.out
#SBATCH --error=logs/script61_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== Script 61: Feature Engineering (leakage-free per-fold gene selection) ==="
echo "Started: $(date)"
Rscript 61_feature_engineering.R 2>&1
echo "Completed: $(date)"

echo "=== Validation ==="
ls -la ../results/staging_classifier/loco_gene_lists/ 2>/dev/null || echo "No loco_gene_lists directory"
echo "Done: $(date)"
