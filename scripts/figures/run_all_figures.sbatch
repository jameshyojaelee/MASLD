#!/bin/bash
#SBATCH --job-name=all_figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/all_figures_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/all_figures_%j.err

# run_all_figures.sbatch — Master orchestrator for all publication figures
# Runs scripts sequentially within each group, capturing per-script pass/fail.
# Usage: sbatch scripts/figures/run_all_figures.sbatch

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
FIG_SCRIPTS="${BASE}/scripts/figures"
LOG_DIR="${BASE}/RNA-seq/logs"
RESULT_DIR="${LOG_DIR}/fig_results_${SLURM_JOB_ID}"

mkdir -p "${LOG_DIR}" "${RESULT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

echo "========================================================"
echo "  MASLD Publication Figures — Full Regeneration"
echo "  Job: ${SLURM_JOB_ID}"
echo "  Start: $(date)"
echo "  Host:  $(hostname)"
echo "========================================================"

# run_fig SCRIPT LABEL
# Runs one figure script; writes PASS/FAIL to result dir.
run_fig() {
  local script="$1"
  local label="$2"
  local log="${LOG_DIR}/fig_${label// /_}_${SLURM_JOB_ID}.log"
  local result_file="${RESULT_DIR}/${label// /_}"

  if [ ! -f "${FIG_SCRIPTS}/${script}" ]; then
    echo "SKIP" > "${result_file}"
    echo "  [SKIP] ${label} (script not found)"
    return
  fi

  echo "  [RUN ] ${label} ..."
  Rscript "${FIG_SCRIPTS}/${script}" > "${log}" 2>&1
  local rc=$?
  if [ "${rc}" -eq 0 ]; then
    echo "PASS" > "${result_file}"
    echo "  [PASS] ${label}"
  else
    echo "FAIL rc=${rc}" > "${result_file}"
    echo "  [FAIL] ${label} (rc=${rc}) — ${log}"
  fi
}

# run_group GROUP_NAME SCRIPT1 LABEL1 SCRIPT2 LABEL2 ...
run_group() {
  local group_name="$1"
  shift
  echo ""
  echo "--- GROUP: ${group_name} ($(date +%H:%M:%S)) ---"
  while [ "$#" -gt 0 ]; do
    local script="$1"
    local label="$2"
    shift 2
    run_fig "${script}" "${label}"
  done
}

# ============================================================
# GROUP 1 (Main Fig 1 — Atlas overview) removed 2026-06-25:
# fig1_atlas_overview.R / fig1_loo_cv_panels.R / fig1_compact.R were retired;
# the combined fig1_atlas_overview.pdf + fig1_compact.pdf are no longer generated.
# Fig 1 is assembled manually in Illustrator from FIG1_DIR/panels/*.pdf.
# ============================================================
# GROUP 2: Main Fig 2-3 — mid-stage (F1-F3) inflection + deconfounding
# ============================================================
run_group "Main Fig2-3 — Switch + Deconfounding" \
  fig2_progression_sex.R               "Fig2 progression_sex" \
  fig2_concordance_atlas.R       "Fig2 concordance_atlas" \
  fig2_disease_progression.R     "Fig2 disease_progression" \
  fig2_compact.R                 "Fig2 compact" \
  fig3_deconfounding.R           "Fig3 deconfounding" \
  fig3_deconvolution.R           "Fig3 deconvolution" \
  fig3_progression.R             "Fig3 progression" \
  fig3_epigenomic_panels.R       "Fig3 epigenomic_panels" \
  fig3_compact.R                 "Fig3 compact"

# ============================================================
# GROUP 3: Main Fig 4 — Causal architecture
# ============================================================
run_group "Main Fig4 — Causal" \
  fig3_regulatory_architecture.R     "Fig4 causal_architecture" \
  fig4_celltype_coloc_heatmap.R  "Fig4 celltype_coloc" \
  fig4_gwas_volcano.R            "Fig4 gwas_volcano" \
  fig4_pharma_panels.R           "Fig4 pharma_panels" \
  fig4_proteomics_panels.R       "Fig4 proteomics_panels" \
  fig4_spatial_panels.R          "Fig4 spatial_panels" \
  fig4_compact.R                 "Fig4 compact"

# ============================================================
# GROUP 4: Main Fig 5 — Convergence matrix (absorbed legacy Fig 6 content)
# ============================================================
# Archived 2026-04-22: 6-fig → 5-fig restructure.
# Legacy fig6_* scripts (therapeutic_windows, pharmacotranscriptomics,
# translation, upset_nonredundancy, convergence_sankey) moved to
# archive/fig6_ditched_2026-04-22/. Their content is folded into the
# fig5_convergence pipeline + figS_therapeutics supplementary family.
run_group "Main Fig5 — Convergence" \
  fig5_sex_dimorphism.R          "Fig5 sex_dimorphism" \
  fig5_translation.R             "Fig5 translation" \
  fig5_convergence.R             "Fig5 evidence_matrix_v2"

# ============================================================
# GROUP 5: S01-S03 — QC, progression, deconvolution
# ============================================================
run_group "Supp S01-S03 — QC/Progression/Deconv" \
  figS01_qc.R                        "S01 qc" \
  figS_positive_controls.R            "S01 positive_controls" \
  figS_combat_sensitivity.R           "S01 combat_sensitivity" \
  fig_cohort_replication.R            "S01 cohort_replication" \
  figS_deg_validation.R               "S02 deg_validation" \
  fig_progression_landscape_consolidated.R "S02 progression_landscape" \
  fig_contrast_landscape.R            "S02 contrast_landscape" \
  fig_c13_fibrosis_confounding.R      "S02 c13_confounding" \
  fig_c13_consolidated.R              "S02 c13_consolidated" \
  fig_within_nash_progression.R       "S02 within_nash" \
  figS_ccc_dynamics.R                 "S02 ccc_dynamics" \
  figS_disease_progression_supp.R     "S02 progression_supp" \
  figS_sex_divergent.R                "S02 sex_divergent" \
  figS_sc_validation.R                "S03 sc_validation"

# ============================================================
# GROUP 6: S04 — COLOC extended
# ============================================================
run_group "Supp S04 — COLOC" \
  figS04_coloc.R                     "S04 coloc_main" \
  figS_coloc_threshold_sensitivity.R "S04 coloc_threshold" \
  figS_coloc_sensitivity.R           "S04 coloc_sensitivity" \
  figS_coloc_method_comparison.R     "S04 coloc_method" \
  figS_twas_coloc_intact.R           "S04 twas_coloc_intact" \
  figS_scdrs.R                       "S04 scdrs" \
  figS_stratified_enrichment.R       "S04 stratified_enrichment" \
  figS_mediation.R                   "S04 mediation" \
  figS_cross_eqtl.R                  "S04 cross_eqtl" \
  figS_susie_comparison.R            "S04 susie_comparison"

# ============================================================
# GROUP 7: S05-S06 — Epigenomic, spatial, cross-species
# ============================================================
run_group "Supp S05-S06 — Epigenomic/Spatial/Cross-species" \
  figS05_gwas_atac.R                 "S05 gwas_atac" \
  figS_spatial.R                     "S05 spatial" \
  figS_spatial_validation.R          "S05 spatial_validation" \
  figS_finemapping_concordance.R     "S05 finemapping_concordance" \
  figS_sex_epigenomic_validation.R   "S05 sex_epigenomic" \
  figS06_cross_species.R             "S06 cross_species" \
  figS_conserved_core.R              "S06 conserved_core"

# ============================================================
# GROUP 8: S07-S08 — ncRNA, subtyping, convergence
# ============================================================
run_group "Supp S07-S08 — ncRNA/Subtyping/Convergence" \
  figS07_ncrna.R                     "S07 ncrna" \
  figS08_nmf_subtyping.R             "S08 nmf_subtyping" \
  figS_convergence.R                 "S08 convergence" \
  figS_evidence_convergence.R        "S08 evidence_convergence" \
  figS_atlas_validation.R            "S08 atlas_validation" \
  figS_depmap_essentiality.R         "S08 depmap" \
  figS_cohort_forest.R               "S08 cohort_forest" \
  figS_functional_activity.R         "S08 functional_activity"

# ============================================================
# GROUP 8b: S09-S10 — Multi-ancestry, prediction, staging
# ============================================================
run_group "Supp S09-S10 — Multi-ancestry/Prediction" \
  figS09_multi_ancestry_coloc.R      "S09 multi_ancestry" \
  figS10_multiprogram.R              "S10 multiprogram" \
  fig_vae_staging.R                  "S10 vae_staging" \
  fig_nas_vae_progression.R          "S10 nas_vae_progression" \
  fig_staging_classifier.R           "S10 staging_classifier" \
  fig_composition_features.R         "S10 composition_features"
  # figS16_plasma_proteomics.R WITHDRAWN 2026-07-04 (round-2 E3c): Olink plasma
  # label-integrity withdrawal; script + inputs archived under
  # data/archive/withdrawn_olink_2026-06-01/. Do not re-add.

# ============================================================
# GROUP 9: Sensitivity figures
# ============================================================
run_group "Supp Sensitivity" \
  figS_sensitivity.R                    "Sens main" \
  figS_published_panel_benchmark.R      "Sens published_benchmark" \
  figS_database_comparison.R            "Sens database_comparison" \
  figS_chen_influence.R                 "Sens chen_influence" \
  figS_bmi_sensitivity.R                "Sens bmi" \
  figS_steatosis_bmi_sensitivity.R      "Sens steatosis_bmi" \
  figS_deg_threshold_landscape.R        "Sens deg_threshold" \
  figS_variance_partition.R             "Sens variance_partition" \
  figS_sex_interaction_subsampling.R    "Sens sex_interaction" \
  figS_cross_disease_specificity.R      "Sens cross_disease"

# ============================================================
# SUMMARY
# ============================================================
echo ""
echo "========================================================"
echo "  SUMMARY  ($(date))"
echo "========================================================"

PASS=0; FAIL=0; SKIP=0
echo "  Per-script results:"
for f in "${RESULT_DIR}"/*; do
  label=$(basename "${f}" | tr '_' ' ')
  result=$(cat "${f}" 2>/dev/null)
  printf "    %-45s %s\n" "${label}" "${result}"
  case "${result}" in
    PASS) PASS=$((PASS+1)) ;;
    SKIP) SKIP=$((SKIP+1)) ;;
    *)    FAIL=$((FAIL+1)) ;;
  esac
done

echo ""
echo "  PASS: ${PASS}  FAIL: ${FAIL}  SKIP: ${SKIP}"

total_pdfs=$(find "${BASE}/figures/main" "${BASE}/figures/supplementary" \
  -name "*.pdf" 2>/dev/null | wc -l)
echo "  Total PDFs in main + supplementary: ${total_pdfs}"
echo ""

if [ "${FAIL}" -gt 0 ]; then
  echo "  WARNING: ${FAIL} script(s) failed. Logs in ${LOG_DIR}/"
  exit 1
fi
echo "  All done!"
