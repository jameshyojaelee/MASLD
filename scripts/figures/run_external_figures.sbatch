#!/bin/bash
#SBATCH --job-name=ext_figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/ext_figures_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/ext_figures_%j.err

# run_external_figures.sbatch
# Regenerates figures from scripts OUTSIDE scripts/figures/ that were not
# covered by run_all_figures.sbatch. Targets:
#   - Stratified causal (RNA-seq/209,210b,212,214,216,218)  → figS04_coloc/
#   - GWAS finemapping figures (src/08)                     → figS05_epigenomic_spatial/
#   - GWAS-ATAC figures (src/57)                            → figS05_epigenomic_spatial/
#   - Proteomics concordance (Analysis/Proteomics)          → figS05_epigenomic_spatial/
#   - Cross-species cell-type concordance (Analysis/CCS)    → figS06_cross_species/
#   - Pseudotime figures (Analysis/SingleCell/307*)         → figS02_progression/

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
LOG_DIR="${BASE}/RNA-seq/logs"
RESULT_DIR="${LOG_DIR}/ext_results_${SLURM_JOB_ID}"
mkdir -p "${LOG_DIR}" "${RESULT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

echo "========================================================"
echo "  External Figure Scripts — Regeneration"
echo "  Job: ${SLURM_JOB_ID}"
echo "  Start: $(date)"
echo "========================================================"

run_script() {
  local script="$1"
  local label="$2"
  local log="${LOG_DIR}/ext_${label// /_}_${SLURM_JOB_ID}.log"
  local result_file="${RESULT_DIR}/${label// /_}"

  if [ ! -f "${BASE}/${script}" ]; then
    echo "SKIP" > "${result_file}"
    echo "  [SKIP] ${label} (not found: ${script})"
    return
  fi

  echo "  [RUN ] ${label} ..."
  Rscript "${BASE}/${script}" > "${log}" 2>&1
  local rc=$?
  if [ "${rc}" -eq 0 ]; then
    echo "PASS" > "${result_file}"
    echo "  [PASS] ${label}"
  else
    echo "FAIL rc=${rc}" > "${result_file}"
    echo "  [FAIL] ${label} (rc=${rc}) — ${log}"
  fi
}

# ============================================================
# GROUP A: Stratified causal → figS04_coloc/
# ============================================================
echo ""
echo "--- GROUP A: Stratified causal figures ($(date +%H:%M:%S)) ---"
run_script "RNA-seq/210b_progression_coloc_plots.R"  "A1 progression_coloc_plots"
run_script "RNA-seq/212_progression_figures.R"       "A2 progression_figures"
run_script "RNA-seq/209_sex_subtype_figures.R"       "A3 sex_subtype_figures"
run_script "RNA-seq/214_spatial_figures.R"           "A4 spatial_coloc_figures"
run_script "RNA-seq/216_pharma_figures.R"            "A5 pharma_figures"
run_script "RNA-seq/218_ashr_lfc_threshold.R"        "A6 ashr_lfc_threshold"

# ============================================================
# GROUP B: GWAS finemapping + ATAC → figS05_epigenomic_spatial/
# ============================================================
echo ""
echo "--- GROUP B: GWAS finemapping + ATAC figures ($(date +%H:%M:%S)) ---"
run_script "GWAS/finemapping/src/08_finemapping_figures.R" "B1 finemapping_figures"
run_script "GWAS/finemapping/src/57_gwas_atac_figures.R"   "B2 gwas_atac_figures"

# ============================================================
# GROUP C: Proteomics → figS05_epigenomic_spatial/
# ============================================================
echo ""
echo "--- GROUP C: Proteomics ($(date +%H:%M:%S)) ---"
run_script "Analysis/Proteomics/scripts/mrna_protein_concordance.R" "C1 mrna_protein_concordance"

# ============================================================
# GROUP D: Cross-species cell-type concordance → figS06_cross_species/
# ============================================================
echo ""
echo "--- GROUP D: Cross-species concordance ($(date +%H:%M:%S)) ---"
run_script "Analysis/Cross_Species_Concordance/scripts/06_celltype_concordance.R" "D1 celltype_concordance"

# ============================================================
# GROUP E: Pseudotime figures → figS02_progression/
# Only run the final consolidated figure (307g); others are intermediate
# ============================================================
echo ""
echo "--- GROUP E: Pseudotime figures ($(date +%H:%M:%S)) ---"
run_script "Analysis/SingleCell/scripts/307g_pseudotime_final.R"  "E1 pseudotime_final"
run_script "Analysis/SingleCell/scripts/307d_heatmap_with_nas.R"  "E2 pseudotime_heatmap"
run_script "Analysis/SingleCell/scripts/307e_virtual_staging.R"   "E3 virtual_staging"

# ============================================================
# SUMMARY
# ============================================================
echo ""
echo "========================================================"
echo "  SUMMARY  ($(date))"
echo "========================================================"

PASS=0; FAIL=0; SKIP=0
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
echo ""

echo "  figS04_coloc PDFs:          $(find "${BASE}/figures/supplementary/figS04_coloc" -name "*.pdf" | wc -l)"
echo "  figS05_epigenomic PDFs:     $(find "${BASE}/figures/supplementary/figS05_epigenomic_spatial" -name "*.pdf" | wc -l)"
echo "  figS06_cross_species PDFs:  $(find "${BASE}/figures/supplementary/figS06_cross_species" -name "*.pdf" | wc -l)"
echo "  figS02_progression PDFs:    $(find "${BASE}/figures/supplementary/figS02_progression" -name "*.pdf" | wc -l)"

if [ "${FAIL}" -gt 0 ]; then
  echo ""
  echo "  WARNING: ${FAIL} script(s) failed. Logs in ${LOG_DIR}/"
  exit 1
fi
echo "  All done!"
