#!/bin/bash
#SBATCH --job-name=regen_figs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/regen_figs_%j.log

set -uo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

SCRIPTS=(
  scripts/figures/fig1_atlas_overview.R
  scripts/figures/fig1_compact.R
  scripts/figures/fig1_integration_value_panels_v2.R
  scripts/figures/fig1_loo_cv_panels.R
  scripts/figures/fig2_compact.R
  scripts/figures/fig2_concordance_atlas.R
  scripts/figures/fig2_disease_progression.R
  scripts/figures/fig3_compact.R
  # fig3_deconv_singlecell.R archived — content promoted to fig2_compact.R
  scripts/figures/fig3_deconvolution.R
  scripts/figures/fig3_epigenomic_panels.R
  scripts/figures/fig4_proteomics_panels.R
  scripts/figures/fig4_spatial_panels.R
  scripts/figures/fig5_causal_architecture.R
  scripts/figures/fig5_translation.R
  # fig6_pharmacotranscriptomics.R — archived 2026-04-22 (6-fig → 5-fig restructure; content folded into fig5_convergence)
  # fig6_translation.R — archived 2026-04-22 (6-fig → 5-fig restructure)
  scripts/figures/fig7_multi_evidence.R
  scripts/figures/fig7_validation.R
  scripts/figures/fig_sc_compositional_crossspecies.R
  scripts/figures/figS_combat_sensitivity.R
  scripts/figures/figS_deg_validation.R
  scripts/figures/figS_positive_controls.R
  scripts/figures/figS_sex_divergent.R
  scripts/figures/figS4_conserved.R
  scripts/figures/figS6_per_cohort_qc.R
  scripts/figures/figS8_disease_progression_supp.R
  scripts/figures/figS_atac.R
  scripts/figures/figS_causal_extended.R
  scripts/figures/figS_convergence.R
  scripts/figures/figS03_deconvolution.R
  scripts/figures/figS_functional_activity.R
  scripts/figures/figS_mouse_overview.R
  scripts/figures/figS_sc_validation.R
  scripts/figures/figS_singlecell.R
  scripts/figures/figS_spatial.R
  scripts/figures/figS_spatial_validation.R
)

TOTAL=${#SCRIPTS[@]}
PASS=0
FAIL=0

for i in "${!SCRIPTS[@]}"; do
  script="${SCRIPTS[$i]}"
  name=$(basename "$script")
  echo "=== [$((i+1))/$TOTAL] $name ==="
  if Rscript "$script" 2>&1; then
    echo "  -> OK"
    ((PASS++))
  else
    echo "  -> FAILED (exit $?)"
    ((FAIL++))
  fi
  echo ""
done

echo "============================================"
echo "DONE: $PASS passed, $FAIL failed out of $TOTAL"
echo "============================================"
