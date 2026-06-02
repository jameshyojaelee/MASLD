#!/usr/bin/env bash
#SBATCH --job-name=compare_methods
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/compare_methods_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/compare_methods_%j.err

# MuSiC vs BayesPrism method comparison across all 8 datasets.
# Generates correlation and bar plot PDFs in comparison_plots/.
# Resources: 32G / 2h for loading all 8 datasets' proportion matrices.
#
# DEPENDENCY: Runs after run_bayesprism_missing_human.sbatch (dependency chain).

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV_PREFIX=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/music_deconv

R_BIN="${ENV_PREFIX}/bin/R"
export R_LIBS_USER=""
export R_LIBS_SITE="${ENV_PREFIX}/lib/R/library"
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

echo "=============================="
echo "Deconvolution methods comparison (all 8 datasets)"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "=============================="

"${R_BIN}" --no-save --no-restore --slave \
  --file="${PROJECT_ROOT}/Analysis/Deconvolution/scripts/16_compare_methods.R"

echo "Finished: $(date)"
echo "Outputs in ${PROJECT_ROOT}/Analysis/Deconvolution/comparison_plots/"
ls -lh "${PROJECT_ROOT}/Analysis/Deconvolution/comparison_plots/" 2>/dev/null || \
  echo "  (directory not created -- check for errors)"
