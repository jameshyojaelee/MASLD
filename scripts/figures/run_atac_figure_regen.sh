#!/bin/bash
#SBATCH --job-name=atac_fig_regen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=scripts/figures/logs/atac_fig_regen_%j.out
#SBATCH --error=scripts/figures/logs/atac_fig_regen_%j.err

# ==========================================================================
# Regenerate all figures affected by scATAC label transfer
#
# Main figures:
#   - Fig 2 panel (f): chromVAR top motifs (now reads chromvar_v2)
#   - Fig 3 panels (f-i): epigenomic panels including chromVAR volcano
#
# Supplementary:
#   - figS_atac Python panels (01-08): UMAP, proportions, markers with new labels
#   - figS_atac R panels (09-16): DA volcano, chromVAR heatmap with v2 data
# ==========================================================================

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$(pwd)"
mkdir -p scripts/figures/logs figures/supplementary/figS05_epigenomic_spatial

echo "============================================================"
echo "Regenerating ATAC-affected figures after label transfer"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Date:   $(date)"
echo "============================================================"

# --- Phase 1: Supplementary Python panels (snapatac2 env) ----------------
echo ""
echo ">>> Phase 1: Supplementary Python panels (01-08)..."
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

python scripts/figures/figS_atac_python_panels.py
echo ">>> Python panels done."

# --- Phase 2: Main + Supplementary R figures (rnaseq env) -----------------
echo ""
echo ">>> Phase 2: R figures (Fig 3 epigenomic, figS_atac)..."
micromamba activate rnaseq

echo ""
echo "--- Fig 3 epigenomic panels ---"
Rscript scripts/figures/fig3_epigenomic_panels.R 2>&1 || echo "WARNING: fig3_epigenomic_panels.R had errors"

echo ""
echo "--- Supplementary ATAC R panels (09-16) ---"
Rscript scripts/figures/figS_atac.R 2>&1 || echo "WARNING: figS_atac.R had errors"

# --- Summary --------------------------------------------------------------
echo ""
echo "============================================================"
echo "Figure regeneration COMPLETE: $(date)"
echo ""
echo "Main figures:"
ls -lh figures/fig2/*.pdf 2>/dev/null || echo "  (no fig2 PDFs found)"
echo ""
echo "Epigenomic panels:"
ls -lh figures/fig3/panels/fig3_panel_*.pdf 2>/dev/null || echo "  (no fig3 epigenomic panels)"
echo ""
echo "Supplementary ATAC:"
ls -lh figures/supplementary/figS05_epigenomic_spatial/*.pdf 2>/dev/null | head -20
echo "============================================================"
