#!/bin/bash
#SBATCH --job-name=liver_pub_figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=scripts/figures/logs/pub_figures_%j.out
#SBATCH --error=scripts/figures/logs/pub_figures_%j.err

set -eo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${BASE}/scripts/figures"
FIG_OUT="${BASE}/figures"
LOG_DIR="${SCRIPT_DIR}/logs"

mkdir -p "${FIG_OUT}" "${LOG_DIR}"

cd "${BASE}"

# Activate environment
eval "$(micromamba shell hook --shell bash)"
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
    echo "[$(date '+%H:%M:%S')] ${name} FAILED"
    failed+=("${name}")
  fi
  echo ""
}

# ---------- Main figures (1-5) ----------
# Fig 1: Multi-Cohort MASLD Transcriptomic Atlas (reverted)
# Fig 2: Single-Cell Atlas & Deconvolution (12 panels a-l)
# Fig 3: Causal Architecture (eQTL, GWAS, COLOC, TWAS, INTACT; MR ditched 2026-04-22)
# Fig 4: Pharmacotranscriptomics, Proteomics & Spatial (sub-scripts)
# Fig 5: Translational Convergence & Validation (was Fig 6; Fig 5 convergence → supplementary)
MAIN_FIGURES=(
  "fig1_compact.R"
  # 2026-05-28 (P0-H): fig2_compact.R was retired in the 2026-05-17 Fig 2
  # redesign (commit 20cf5af archived it to
  # archive/fig2_pre_redesign_2026-05-17/fig2_compact.R). Fig 2 is now built
  # from ~22 fig2_panel_*.R scripts that each render a panel PDF into
  # FIG2_DIR/panels/, then tiled by fig2_progression.R (driven by the
  # run_assemble_fig2.sh wrapper). The single-script compositor no longer
  # exists, so fig2_progression.R is the correct compositor to call here.
  "fig2_progression.R"
  "fig3_compact.R"
  # fig4_compact.R sources pharma + spatial sub-scripts and runs the canonical
  # fig4_validation.R first (proteomics + spatial composite).
  "fig4_compact.R"
  "fig5_translation.R"
)

for script_name in "${MAIN_FIGURES[@]}"; do
  script="${SCRIPT_DIR}/${script_name}"
  if [[ -f "${script}" ]]; then
    run_figure "${script}"
  else
    echo "WARNING: ${script_name} not found"
  fi
done

# ---------- Supplementary figures (7) ----------
# S1: Sunburst + UpSet + cumulative discovery curve
# S2: ComBat-seq sensitivity
# S3: Sex inference validation + sex-divergent sensitivity
# S4: BayesPrism cross-method + pseudobulk-bulk concordance
# S5: SCENIC+/chromVAR extended + extended spatial
# S6: Multi-ancestry COLOC extended
# S7: Positive controls + variance partition + per-cohort QC
SUPP_FIGURES=(
  "figS_convergence.R"
  "figS_atlas_validation.R"
  "figS_combat_sensitivity.R"
  "figS_positive_controls.R"
  "figS_sex_divergent.R"
  # 2026-05-28 (P0-H): stale names fixed. The conserved-enrichment supplement
  # is figS_conserved_core.R (was referenced as figS4_conserved.R) and the
  # variance-partition supplement is figS_variance_partition.R (was referenced
  # as figS5_variance_partition.R). Both real scripts exist; only the runner
  # references were out of date.
  "figS_conserved_core.R"
  "figS_variance_partition.R"
  "figS6_per_cohort_qc.R"
  "supp_fig7_extended_celltype.R"
  "figS_cohort_forest.R"
  "figS_spatial.R"
  "figS_spatial_validation.R"
  # 2026-05-28 (P0-H): figS_causal_extended.R removed (deleted 2026-04-09,
  # superseded). Its three panels are now covered by scripts already in this
  # list: (a) cross-ancestry PP.H4 scatter -> figS_cross_eqtl.R + figS04_coloc.R;
  # (b) COLOC prior-sensitivity heatmap -> figS_coloc_sensitivity.R;
  # (c) COLOC GWAS-overlap UpSet -> figS04_coloc.R + figS_susie_comparison.R.
  "figS_sensitivity.R"
  "figS_coloc_sensitivity.R"
  "figS_susie_comparison.R"
  "figS_cross_eqtl.R"
)

for script_name in "${SUPP_FIGURES[@]}"; do
  script="${SCRIPT_DIR}/${script_name}"
  if [[ -f "${script}" ]]; then
    run_figure "${script}"
  else
    echo "WARNING: ${script_name} not found"
  fi
done

# ---------- Summary ----------
echo "=========================================="
echo "FIGURE GENERATION COMPLETE"
echo "=========================================="
echo ""
echo "Output directory: ${FIG_OUT}"
echo ""
echo "Generated files:"
ls -lh "${FIG_OUT}"/fig*.pdf 2>/dev/null || echo "  (no PDFs found)"
echo ""
ls -lh "${FIG_OUT}"/fig*.png 2>/dev/null || echo "  (no PNGs found)"
echo ""

if [[ ${#failed[@]} -gt 0 ]]; then
  echo "FAILED scripts (${#failed[@]}):"
  for f in "${failed[@]}"; do
    echo "  - ${f}"
  done
  exit 1
else
  echo "All figures generated successfully!"
  exit 0
fi
