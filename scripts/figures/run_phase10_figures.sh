#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=phase10_figs
#SBATCH --output=scripts/figures/logs/phase10_figs_%A_%a.out
#SBATCH --error=scripts/figures/logs/phase10_figs_%A_%a.err
#SBATCH --array=0-22

# Phase 10 Team C: re-render all main + supplementary figures against the
# PolyFun-backed canonical SuSiE-COLOC atlas.
#
# The fig3 panel script (index 1) must complete before assemble_fig3.R runs;
# assemble_fig3.R is invoked separately AFTER the array completes.

set -eo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${PROJECT_ROOT}"

# Activate environment (set -u disabled because conda/micromamba activate
# scripts reference unbound vars during init).
source ~/.bashrc
micromamba activate rnaseq

# Figure script registry (index = SLURM_ARRAY_TASK_ID)
SCRIPTS=(
  # Main figures (0-2)
  "scripts/figures/fig1_atlas_overview_v2.R"
  "scripts/figures/fig3_regulatory_architecture_v2.R"
  "scripts/figures/fig5_convergence_v3.R"
  # Supplementary COLOC (4-8)
  "scripts/figures/figS04_coloc.R"
  "scripts/figures/figS05_gwas_atac.R"
  "scripts/figures/figS09_multi_ancestry_coloc.R"
  "scripts/figures/figS09_locus_zoom.R"
  "scripts/figures/figS10_multiprogram.R"
  # Convergence + sensitivity (9-14)
  "scripts/figures/figS_convergence.R"
  "scripts/figures/figS_convergence_evidence.R"
  "scripts/figures/figS_coloc_sensitivity.R"
  "scripts/figures/figS_coloc_threshold_sensitivity.R"
  "scripts/figures/figS_coloc_method_comparison.R"
  "scripts/figures/figS_coloc_window_sensitivity.R"
  # SuSiE validation (15-18)
  "scripts/figures/figS_susie_expression_scatters.R"
  "scripts/figures/figS_susie_gwas_volcano.R"
  "scripts/figures/figS_susie_celltype_coloc_heatmap.R"
  "scripts/figures/figS_susie_comparison.R"
  # Stratified plots (RNA-seq/) (19-23)
  "RNA-seq/209_sex_subtype_figures.R"
  "RNA-seq/210b_progression_coloc_plots.R"
  "RNA-seq/212_progression_figures.R"
  "RNA-seq/214_spatial_figures.R"
  "RNA-seq/216_pharma_figures.R"
)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
SCRIPT="${SCRIPTS[$IDX]}"

if [[ -z "${SCRIPT:-}" ]]; then
  echo "[ERROR] No script registered for array index ${IDX}" >&2
  exit 1
fi

echo "============================================================"
echo "Phase 10 Team C: Re-rendering ${SCRIPT}"
echo "Array index: ${IDX}"
echo "Job ID: ${SLURM_JOB_ID:-N/A}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

if [[ ! -f "${SCRIPT}" ]]; then
  echo "[ERROR] Script not found: ${SCRIPT}" >&2
  exit 2
fi

# Run from project root so relative path constants in load_figure_data.R resolve.
Rscript "${SCRIPT}"
RC=$?

echo "============================================================"
echo "Finished: $(date)"
echo "Exit code: ${RC}"
echo "============================================================"

exit ${RC}
