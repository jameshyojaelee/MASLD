#!/bin/bash
#SBATCH --job-name=figregen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/rerender_gse_labels_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/rerender_gse_labels_%j.err
#
# Re-render every figure script whose author-name display labels were swapped to
# GEO/accession IDs (author->GSE sweep). Continue-on-failure with per-script
# PASS/FAIL + duration logging so one broken/data-missing script does not block
# the rest. Render failures here are environmental/data (edits were display-label
# only and all parse-clean), and are reported, not silently swallowed.

# NB: do NOT use `set -u` — micromamba/conda activate.d scripts reference unbound
# vars (e.g. ADDR2LINE) and abort the whole job under nounset.
set -o pipefail
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SF="${BASE}/scripts/figures"
export MASLD_PROJECT_ROOT="${BASE}"
cd "${BASE}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

PER_SCRIPT_TIMEOUT=1800   # 30 min hard cap per script
PASS=(); FAIL=()

run() {  # run <interpreter> <relative-script-and-args...>
  local kind="$1"; shift
  local disp="$*"
  echo ""
  echo "===================================================================="
  echo ">>> [$kind] $disp"
  echo "===================================================================="
  local t0=$SECONDS
  if [ "$kind" = "R" ]; then
    timeout "${PER_SCRIPT_TIMEOUT}" Rscript "${SF}/$1"
  else
    timeout "${PER_SCRIPT_TIMEOUT}" python "${SF}/$@"
  fi
  local rc=$?
  local dt=$((SECONDS - t0))
  if [ $rc -eq 0 ]; then PASS+=("$disp (${dt}s)"); echo "[PASS] $disp (${dt}s)";
  else FAIL+=("$disp (rc=$rc, ${dt}s)"); echo "[FAIL] $disp (rc=$rc, ${dt}s)"; fi
}

############################ Python overview figures ############################
run PY fig_bulkrna_matrix.py --panel
run PY metadata_matrix.py
run PY dataset_treemap.py
run PY dataset_overview_figure.py
run PY published_deg_comparison.py

################################ R figures ######################################
R_SCRIPTS=(
  fig1_umap.R
  fig1_atlas_overview_v2.R
  figS_batch_harmony_sweep.R
  figS_batch_correction.R
  figS_batch_correction_extra.R
  figS_batch_correction_overlap.R
  figS01_batch_correction_pca.R
  figS_pca_definitive.R
  figS_pca_clinical.R
  figS01_mitochondrial_pct.R
  figS_multimethod_batch_pca.R
  figS_multimethod_batch_pca_split.R
  figS_multimethod_batch_pca_split_top200.R
  figS_multimethod_disease_separation.R
  figS_pls_pooled.R
  figS_pls_disease_separation.R
  figS_supervised_disease_pca_top100DEG.R
  figS_supervised_disease_pca_top200DEG.R
  figS_supervised_disease_pca_top200DEG_clinical.R
  figS_disease_signature_sweep.R
  figS_chen_influence.R
  figS_integration_value.R
  figS_integration_value_venn.R
  loo_cv_stability.R
  loo_cv_effect_preservation.R
  figS_robust_loco.R
  figS_lfc_sensitivity_loo.R
  figS_lfc_sensitivity_loo_lines.R
  figS_lfc_sensitivity_loo_robustness.R
  figS_lfc_sensitivity_loo_mash_masl.R
  figS_lfc_sensitivity_loo_lfc_compare.R
  figS_lfc_sensitivity_loo_concordance.R
  per_study_lfc_correlation_heatmap.R
  figS_sex_dimorphism.R
  deg_upset_cohorts.R
  fig4_spatial_panels.R
  cyp3a4_zonation.R
  fig5_panel_calibration.R
  spatial_svg_dynamics.R
  figS_spatial.R
  figS_spatial_validation.R
)
for s in "${R_SCRIPTS[@]}"; do run R "$s"; done

############################### Summary #########################################
echo ""
echo "####################################################################"
echo "RE-RENDER SUMMARY:  ${#PASS[@]} PASS / ${#FAIL[@]} FAIL"
echo "####################################################################"
echo "--- PASS ---"; printf '  %s\n' "${PASS[@]:-（none）}"
echo "--- FAIL ---"; printf '  %s\n' "${FAIL[@]:-（none）}"
