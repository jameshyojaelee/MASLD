#!/bin/bash
#SBATCH --job-name=regen_deg_figs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/regen_deg_figs_%j.out
#SBATCH --error=scripts/figures/logs/regen_deg_figs_%j.err

# Regenerate supplementary + fig6 figures affected by DEG definition change:
#   padj<0.05 + |logFC|>0.3 (5,484 DEGs)

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${BASE}/scripts/figures"
LOG_DIR="${SCRIPT_DIR}/logs"

mkdir -p "${LOG_DIR}"
cd "${BASE}"

# Activate environment
source ~/.bashrc
micromamba activate rnaseq

# Install ggrastr if missing (for rasterized scatter plots)
Rscript -e 'if (!requireNamespace("ggrastr", quietly=TRUE)) install.packages("ggrastr", repos="https://cloud.r-project.org", quiet=TRUE)' 2>/dev/null || true

export MASLD_PROJECT_ROOT="${BASE}"

failed=()

run_figure() {
  local script="$1"
  local name=$(basename "$script" .R)
  echo "=========================================="
  echo "[$(date '+%H:%M:%S')] Running ${name}..."
  echo "=========================================="
  if Rscript "${script}"; then
    echo "[$(date '+%H:%M:%S')] ${name} SUCCEEDED"
  else
    echo "[$(date '+%H:%M:%S')] ${name} FAILED (exit $?)"
    failed+=("${name}")
  fi
  echo ""
}

# ---- Supplementary figures that use dream DEG thresholds ----
SCRIPTS=(
  # Core supplementary figures (all modified for new DEG definition)
  "figS_sensitivity.R"
  "figS_depmap_essentiality.R"
  "figS_causal_extended.R"
  "figS_ncrna.R"
  "figS_evidence_convergence.R"
  "figS_convergence.R"
  "figS_atlas_validation.R"
  "figS_sc_validation.R"
  "figS_sex_subsampling.R"
  "figS_cohort_forest.R"
  "figS_chen_influence.R"
  # Additional supplementary (load dream but use own thresholds — re-run for consistency)
  "figS_mouse_overview.R"
  "figS_atac.R"
  "figS_progression_landscape.R"
  "figS_stage_transitions.R"
  # Main figures (fig5 uses atlas/dream DEG thresholds)
  "fig5_causal_architecture.R"
  "fig5_translation.R"
  # fig6_therapeutic_windows.R and fig6_convergence_sankey.R archived 2026-04-22
  # (6-fig → 5-fig restructure; content folded into fig5_convergence.R + figS_therapeutics)
  "fig5_convergence.R"
  # Fig1-4: use canonical bulk logFC (raw) via load_figure_data.R — need refresh after column rename
  "fig1_integration_value_panels_v2.R"
  "fig2_compact.R"
  "fig3_compact.R"
  "fig3_progression.R"
  "fig4_compact.R"
  "fig4_celltype_coloc_heatmap.R"
)

for script_name in "${SCRIPTS[@]}"; do
  script="${SCRIPT_DIR}/${script_name}"
  if [[ -f "${script}" ]]; then
    run_figure "${script}"
  else
    echo "WARNING: ${script_name} not found, skipping"
  fi
done

# ---- Summary ----
echo "=========================================="
echo "DEG FIGURE REGENERATION COMPLETE"
echo "=========================================="
echo ""
echo "Figures output to: ${BASE}/figures/"
echo ""

if [[ ${#failed[@]} -gt 0 ]]; then
  echo "FAILED scripts (${#failed[@]}):"
  for f in "${failed[@]}"; do
    echo "  - ${f}"
  done
  exit 1
else
  echo "All ${#SCRIPTS[@]} figures regenerated successfully!"
  exit 0
fi
