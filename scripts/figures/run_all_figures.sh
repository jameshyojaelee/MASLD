#!/bin/bash
#SBATCH --job-name=all_pub_figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=scripts/figures/logs/all_figures_%j.out
#SBATCH --error=scripts/figures/logs/all_figures_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${BASE}/scripts/figures"
FIG_OUT="${BASE}/figures"
LOG_DIR="${SCRIPT_DIR}/logs"

mkdir -p "${FIG_OUT}" "${LOG_DIR}"

cd "${BASE}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

# Install ggrastr if missing (for rasterized scatter plots)
Rscript -e 'if (!requireNamespace("ggrastr", quietly=TRUE)) install.packages("ggrastr", repos="https://cloud.r-project.org", quiet=TRUE)' 2>/dev/null || true

export MASLD_PROJECT_ROOT="${BASE}"

failed=()
succeeded=()

run_figure() {
  local script="$1"
  local name=$(basename "$script" .R)
  echo "=========================================="
  echo "[$(date '+%H:%M:%S')] Running ${name}..."
  echo "=========================================="
  if Rscript "${script}" 2>&1; then
    echo "[$(date '+%H:%M:%S')] ${name} SUCCEEDED"
    succeeded+=("${name}")
  else
    echo "[$(date '+%H:%M:%S')] ${name} FAILED"
    failed+=("${name}")
  fi
  echo ""
}

# ---------- All figure scripts ----------
# Skip fig3_deconvolution.R (already generated, needs h5ad)
# Skip supp_fig7_extended_celltype.R (already generated)

FIGURE_SCRIPTS=(
  # Main figures (consolidated)
  "fig2_concordance_atlas.R"
  "fig3_deconvolution.R"
  # fig4_sex_stratification.R — panels moved into fig2 (i, j)
  "fig5_causal_architecture.R"
  "figS_singlecell.R"
  # fig6_pharmacotranscriptomics.R — archived 2026-04-22 (6-fig → 5-fig restructure; content absorbed into fig5_convergence)
  "fig7_multi_evidence.R"
  # Supplementary figures
  "figS_deg_validation.R"
  "figS_combat_sensitivity.R"
  "figS_positive_controls.R"
  "figS_sex_divergent.R"
  "figS4_conserved.R"
  "figS5_variance_partition.R"
  "figS6_per_cohort_qc.R"
  "supp_fig7_extended_celltype.R"
  "fig2_disease_progression.R"
  "fig7_validation.R"
  "figS_cohort_forest.R"
  "figS_mouse_overview.R"
  "figS_evidence_convergence.R"
  "figS_sc_validation.R"
  "figS_spatial.R"
)

for script_name in "${FIGURE_SCRIPTS[@]}"; do
  script="${SCRIPT_DIR}/${script_name}"
  if [[ -f "${script}" ]]; then
    run_figure "${script}"
  else
    echo "WARNING: ${script_name} not found, skipping"
  fi
done

# ---------- Summary ----------
echo "=========================================="
echo "FIGURE GENERATION COMPLETE"
echo "=========================================="
echo ""
echo "Output directory: ${FIG_OUT}"
echo ""
echo "SUCCEEDED (${#succeeded[@]}):"
for s in "${succeeded[@]}"; do
  echo "  + ${s}"
done
echo ""
echo "Generated PDFs:"
ls -lh "${FIG_OUT}"/*.pdf 2>/dev/null || echo "  (no PDFs found)"
echo ""
echo "Generated PNGs:"
ls -lh "${FIG_OUT}"/*.png 2>/dev/null || echo "  (no PNGs found)"
echo ""

if [[ ${#failed[@]} -gt 0 ]]; then
  echo "FAILED (${#failed[@]}):"
  for f in "${failed[@]}"; do
    echo "  - ${f}"
  done
  exit 1
else
  echo "All figures generated successfully!"
  exit 0
fi
