#!/bin/bash
#SBATCH --job-name=figs_final_17gwas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/figs_final_17gwas_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/figs_final_17gwas_%j.err

set -euo pipefail

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${MASLD_PROJECT_ROOT}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

run_fig() {
  local script="$1"
  local label="$2"
  echo "--- ${label} ($(date +%H:%M:%S)) ---"
  Rscript "scripts/figures/${script}" && echo "  OK: ${label}" || echo "  FAIL: ${label}"
  echo ""
}

echo "=== Final figures: 17 EUR GWAS (FinnGen restored) ==="
echo "Start: $(date)"

run_fig figS_coloc_threshold_sensitivity.R "Threshold sensitivity"
run_fig figS04_coloc.R "COLOC ABF overview (figS_coloc_abf.pdf)"
run_fig figS_twas_coloc_intact.R "TWAS x COLOC x INTACT"
run_fig figS_coloc_method_comparison.R "Method comparison"
# RETIRED 2026-07-07: fig4_expression_pip_scatter.R dropped (stale panels, script early-quits)
run_fig figS_coloc_sensitivity.R "Reviewer defense"
run_fig figS_sensitivity.R "Threshold sweep"
run_fig figS_finemapping_concordance.R "Finemapping concordance"

echo ""
echo "=== Atlas rebuild ==="
Rscript RNA-seq/27a_assemble_evidence_atlas.R
Rscript RNA-seq/27b_benchmark_presets.R

echo ""
echo "=== All done: $(date) ==="
