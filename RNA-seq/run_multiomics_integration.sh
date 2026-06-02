#!/bin/bash
#SBATCH --job-name=multi_omic_integration
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/multiomics_%j.out
#SBATCH --error=logs/multiomics_%j.err

# Multi-omic integration analysis: decoupleR + COSMOS
# Runs scripts 51 (decoupleR functional activity) and 52 (COSMOS mechanistic paths)

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p logs

echo "=== Starting multi-omic integration analysis ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-local}"

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

# Step 1: decoupleR functional activity space
echo ""
echo "=== Step 1: decoupleR functional activity ==="
Rscript 51_decoupler_functional_activity.R 2>&1 | tee logs/51_decoupler_${SLURM_JOB_ID:-local}.log

# Step 2: COSMOS mechanistic paths
echo ""
echo "=== Step 2: COSMOS mechanistic paths ==="
Rscript 52_cosmos_mechanistic_paths.R 2>&1 | tee logs/52_cosmos_${SLURM_JOB_ID:-local}.log

echo ""
echo "=== Multi-omic integration complete ==="
echo "Date: $(date)"
echo "Results in: RNA-seq/results/multi_evidence/functional_activity/"
echo "            RNA-seq/results/multi_evidence/cosmos_mechanistic/"
