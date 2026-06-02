#!/bin/bash
#SBATCH --job-name=liver_gwas_overlay
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/gwas_overlay_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/gwas_overlay_%j.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00

set -euo pipefail

# Activate environment
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/figures"

echo "=== Strategy 12: GWAS Spatial Overlay ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-4}"
echo ""

# Ensure output directories exist
mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/gwas_spatial_convergence
mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs
mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/figures

Rscript "$SCRIPT_DIR/21_gwas_spatial_overlay.R"

echo ""
echo "=== Done: $(date) ==="
