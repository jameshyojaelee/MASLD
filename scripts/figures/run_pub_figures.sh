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
  # Fig 2/3 (RNA-seq) is assembled manually from the panel PDFs in
  # FIG2_DIR/panels/ — no compositor script (removed 2026-06-11).
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

# ---------- Fig 3 (RNA-seq) individual panels ----------
# Fig 3 / S3 is assembled from individual panel scripts (no compositor). Added
# 2026-06-17 (main 3H-3M + supp S3J-S3T); updated 2026-06-18 (+ 3N single-cell
# cascade; zonation 3J/S3R relocated to Fig 4, see FIG4_SPATIAL_PANELS below;
# #2 fibrosis mechanism cascade deleted). Callout->filename index in
# figures/main/fig3_RNAseq/README.md. Plotting-only (main -> fig3_RNAseq/panels/,
# supp -> figures/supplementary/<theme dir>/).
FIG3_RNASEQ_PANELS=(
  # main Fig 3 (-> figures/main/fig3_RNAseq/panels/)
  "progression_cascade.R"                              # 3E (comprehensive stage-progression cascade)
  "drug_target_orthogonality_landscape.R"              # 3H
  "crossmodal_convergence_matrix.R"                    # 3I (spatial col dropped 2026-06-18)
  "nmf_dominant_program_stage_composition.R"           # 3K
  "multicelltype_hotspot_cascade_landscape.R"          # 3M
  "singlecell_disease_cascade.R"                       # 3N (shared-x single-cell cascade)
  # supplementary (-> figures/supplementary/<theme dir>/)
  "per_study_lfc_correlation_heatmap.R"                # S3J
  "deg_meta_heterogeneity_ridge.R"                     # S3K
  "loo_per_gene_reproducibility_bar.R"                 # S3L
  "fibrosis_stage_directional_asymmetry.R"             # S3M
  "fibrosis_stage_hallmark_escalation.R"               # S3N
  "deconvolution_composition_shift_forest.R"           # S3O
  "scrna_clr_abundance_forest.R"                       # S3P
  "hep_module_sc_vs_bulk_concordance_scatter.R"        # S3Q
  "progression_driver_coloc_phenotype_class_heatmap.R" # S3S
  "cross_species_pathway_translatability_matrix.R"     # S3T
  "figS_pca_definitive.R"                              # S3U (definitive control-vs-disease PCA; loads merged DGE -> heavier)
)

for script_name in "${FIG3_RNASEQ_PANELS[@]}"; do
  script="${SCRIPT_DIR}/${script_name}"
  if [[ -f "${script}" ]]; then
    run_figure "${script}"
  else
    echo "WARNING: ${script_name} not found"
  fi
done

# ---------- Fig 4 (spatial/zonation) panels relocated from fig3 (2026-06-18) ----------
# Render to figures/main/fig4_validation/ (scripts repointed to FIG4_DIR); moved
# out of fig3 under the "all spatial -> Fig 4" rule.
FIG4_SPATIAL_PANELS=(
  "zonation_directional_polarity.R"
  "zonation_directional_polarity_xspecies.R"
  "zonation_crossspecies_concordance.R"
)
for script_name in "${FIG4_SPATIAL_PANELS[@]}"; do
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
