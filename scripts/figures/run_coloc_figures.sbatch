#!/bin/bash
#SBATCH --job-name=coloc_figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/coloc_figures_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/coloc_figures_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FIGS="${BASE}/scripts/figures"
export MASLD_PROJECT_ROOT="${BASE}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== COLOC Figure Regeneration ==="
echo "Job: ${SLURM_JOB_ID} | Node: $(hostname) | Start: $(date)"
echo ""

run_fig() {
  local script="$1"
  local label="$2"
  echo "--- ${label} ($(date +%H:%M:%S)) ---"
  Rscript "${FIGS}/${script}" && echo "  OK: ${label}" || echo "  FAILED: ${label}"
  echo ""
}

# New threshold sensitivity figure (from audit)
run_fig "figS_coloc_threshold_sensitivity.R" "figS_coloc_threshold_sensitivity"

# Updated ABF COLOC figure (MHC/LD cluster annotations)
run_fig "figS_coloc_abf.R" "figS_coloc_abf"

# TWAS + COLOC + INTACT figure
run_fig "figS_twas_coloc_intact.R" "figS_twas_coloc_intact"

# Method comparison (ABF vs SuSiE concordance)
run_fig "figS_coloc_method_comparison.R" "figS_coloc_method_comparison"

# RETIRED 2026-07-07: fig4_expression_pip_scatter.R dropped — its outputs
# (panel_expression_pip_scatter/panel_expression_coloc_scatter) are stale panels
# not in the Fig2/FigS2 set; the script itself now early-quits.

# Existing reviewer defense sensitivity figure
run_fig "figS_coloc_sensitivity.R" "figS_coloc_sensitivity"

echo "=== Done: $(date) ==="
