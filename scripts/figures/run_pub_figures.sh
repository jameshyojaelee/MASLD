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
  # Fig 1 (atlas overview) is assembled manually in Illustrator from the panel
  # PDFs in FIG1_DIR/panels/ — no compositor script (fig1_compact.R removed
  # 2026-06-25; the combined fig1_atlas_overview.pdf / fig1_compact.pdf are no
  # longer generated).
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
  # figS_conserved_core.R REMOVED 2026-07-05 — cross-species conserved-core dropped
  # from the atlas paper (no cross-species); script moved to scripts/figures/_legacy/.
  # variance-partition supplement is figS_variance_partition.R (real script exists).
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
  # A–J layout 2026-07-02. Assembled manually in Illustrator (no compositor). Reading order:
  # row1 = 3A metadata / 3B fib×NAS grid / 3C alluvial; then 3D PCA · 3E disconnect ·
  # 3F cascade · 3G UMAP · 3H TF (aligns under 3F, shared F0–F4 axis) · 3I module · 3J L-R.
  # 3A cohort metadata = fig_bulkrna_matrix.py --panel (PYTHON; run via run_regen_fig1_fig5.sh, NOT this R array).
  "nas_fib_grid.R"                                     # 3B fibrosis × NAS patient-count grid
  "figS_integration_value.R"                           # 3C DEG-overlap alluvial (also writes the supp source copy)
  # 3D = figS_pca_definitive.R — listed once in the supplementary block below (emits BOTH fig3d + S3U)
  "coloc_deg_bridge.R"                                 # 3E genetics<->expression disconnect (gene-level bulk LFC vs COLOC PP.H4; promoted from supp 2026-07-02; the 7-lineage cell-type scatter was retired — non-sig correlation)
  "progression_cascade.R"                              # 3F fibrosis-only F0-F4 cascade (shares x-axis with 3H)
  "gen_scrna_umap_embeddable.R"                        # 3G scRNA atlas UMAP
  "tf_convergence_vsF0.R"                              # 3H TF cross-modal convergence (aligns in 3F's column)
  "singlecell_module_heatmap.R"                        # 3I module heatmap + 3J LIANA L-R heatmap (one script)
  "stagedeg_carrier_routing.R"                         # supp figs3_stagedeg_carrier_routing (demoted from main 2026-07-02)
  "nas_activity_cascade.R"                             # supp figs3_nas_activity_cascade (NAS analogue of 3F, in-dir)
  # retired / stale array entries removed 2026-07-02:
  # "genetics_expression_disconnect.R"                 # cell-type disconnect scatter RETIRED (-> _legacy; 7-lineage correlation non-sig p=0.35, hepatocyte-outlier-driven). 3E is now the gene-level coloc_deg_bridge.
  # "deconvolution_celltype_de_limitation.R"           # RETIRED (-> _legacy) with the cell-type scatter it supported. Coarse recompute scripts kept as exploratory: Analysis/SingleCell/scripts/pseudobulk_de_coarse_lineage.R, RNA-seq/201b_celltype_heritability_coarse.R, RNA-seq/89b_toast_celltype_de_coarse.R.
  # "combined_nas_fibrosis_cascade.R"                  # DEPRECATED to _legacy 2026-07-08 — progression_cascade.R now builds the combined NAS+fibrosis cascade (Fig 3E)
  # "drug_target_orthogonality_landscape.R"            # stale entry (not a current fig3 panel)
  # "crossmodal_convergence_matrix.R"                  # stale entry (not a current fig3 panel)
  # "nmf_dominant_program_stage_composition.R"         # (retired)
  # "multicelltype_hotspot_cascade_landscape.R"        # (retired 2026-07-02)
  # "singlecell_disease_cascade.R"                     # (retired 2026-07-02 — line cascade removed)
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
  # cross_species_pathway_translatability_matrix.R (S3T) REMOVED 2026-07-05 — no cross-species; -> _legacy/
  "figS_pca_definitive.R"                              # 3D (fig3d_pca_fibrosis_gradient) + S3U supp variants (loads merged DGE -> heavier)
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
# Render to figures/main/fig5_molecular_context/ (scripts repointed to FIG4_DIR); moved
# out of fig3 under the "all spatial -> Fig 4" rule.
FIG4_SPATIAL_PANELS=(
  "zonation_directional_polarity.R"
  # zonation_directional_polarity_xspecies.R + zonation_crossspecies_concordance.R
  # REMOVED 2026-07-05 — cross-species dropped from the atlas paper; scripts -> _legacy/.
)
for script_name in "${FIG4_SPATIAL_PANELS[@]}"; do
  script="${SCRIPT_DIR}/${script_name}"
  if [[ -f "${script}" ]]; then
    run_figure "${script}"
  else
    echo "WARNING: ${script_name} not found"
  fi
done

# ---------- Fig 4 (validation) — mRNA-protein composite, panel 4c (2026-07-05) ----------
# The one composite kept after the adversarial review; the other five are in
# scripts/figures/_legacy/ (see _legacy/README_composites_cut_2026-07-05.md). Emits to FIG4_DIR.
FIG4_COMPOSITE_PANELS=(
  "composite_mrna_protein.R"   # 4c: protein co-abundance (Reactome-grouped) + mRNA-protein buffering
)
for script_name in "${FIG4_COMPOSITE_PANELS[@]}"; do
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
