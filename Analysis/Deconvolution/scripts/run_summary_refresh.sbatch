#!/usr/bin/env bash
#SBATCH --job-name=deconv_summary
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/summary_refresh_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/summary_refresh_%j.err

# Cross-dataset summary refresh (Phase 3).
# Regenerates creative cross-dataset plots and summary report for all 8 datasets.
# Resources: 64G / 3h provides margin for loading all datasets + plot generation.
#
# Runs:
#   1. 10_make_creative_plots.R  -> results/summary_plots/
#   2. 18_generate_summary.R     -> results/summary_report.md
#
# DEPENDENCY: Runs after QC + method comparison complete (dependency chain).

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV_PREFIX=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/music_deconv

RESULTS_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/results"
BULK_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/bulk"
SCRIPTS_DIR="${PROJECT_ROOT}/Analysis/Deconvolution/scripts"
SUMMARY_REPORT="${RESULTS_BASE}/summary_report.md"

R_BIN="${ENV_PREFIX}/bin/R"
export R_LIBS_USER=""
export R_LIBS_SITE="${ENV_PREFIX}/lib/R/library"
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

echo "=============================="
echo "Cross-dataset summary refresh"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "=============================="

echo ""
echo "--- Step 1: 10_make_creative_plots.R ---"
"${R_BIN}" --no-save --no-restore --slave \
  --file="${SCRIPTS_DIR}/10_make_creative_plots.R" \
  --args \
  "${RESULTS_BASE}" \
  "${BULK_BASE}" \
  "music"

echo ""
echo "--- Step 2: 18_generate_summary.R ---"
"${R_BIN}" --no-save --no-restore --slave \
  --file="${SCRIPTS_DIR}/18_generate_summary.R" \
  --args \
  "${RESULTS_BASE}" \
  "${SUMMARY_REPORT}"

echo ""
echo "Finished: $(date)"
echo "Summary plots in ${RESULTS_BASE}/summary_plots/"
ls -lh "${RESULTS_BASE}/summary_plots/" 2>/dev/null || echo "  (summary_plots/ not created)"
echo "Summary report: ${SUMMARY_REPORT}"
