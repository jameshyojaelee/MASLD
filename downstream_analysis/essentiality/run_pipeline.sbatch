#!/bin/bash
#SBATCH --job-name=essentiality_pipeline
#SBATCH --output=downstream_analysis/essentiality/logs/pipeline_%j.out
#SBATCH --error=downstream_analysis/essentiality/logs/pipeline_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu

# Master Essentiality Analysis Pipeline
# ======================================
# Usage:
#   sbatch downstream_analysis/essentiality/run_pipeline.sbatch

set -euo pipefail

# Environment setup
MICROMAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
ENV_PREFIX=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq

# Enable user-site packages (needed for pip install --user packages like openpyxl)
export PYTHONNOUSERSITE=0
export MPLBACKEND=Agg

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/essentiality

# Get DEG file from argument or use default (relative to project root if needed, or absolute)
# The scripts expect the DEG file path. run_essentiality_pipeline.py defaults to ../../../final_core_degs.csv if not provided?
# The scripts seem to use hardcoded paths in their headers (CORE_DEGS_FILE = ...).
# Let's run them without arguments, assuming they find the files (we verified they do).

echo "========================================"
echo "Essentiality Analysis Pipeline"
echo "========================================"
echo "Start time: $(date)"
echo "Environment: $ENV_PREFIX"
echo "========================================"

echo "[1/3] Running Main Pipeline..."
"${MICROMAMBA}" run -p "${ENV_PREFIX}" python3 scripts/run_essentiality_pipeline.py

echo "[2/3] Running Extended Analysis..."
"${MICROMAMBA}" run -p "${ENV_PREFIX}" python3 scripts/extended_essentiality_analysis.py

echo "[3/3] Running Detailed Analysis (Density/Raincloud)..."
"${MICROMAMBA}" run -p "${ENV_PREFIX}" python3 scripts/analyze_essentiality.py

echo "========================================"
echo "End time: $(date)"
echo "========================================"
