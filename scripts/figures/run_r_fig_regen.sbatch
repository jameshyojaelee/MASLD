#!/bin/bash
#SBATCH --job-name=r_fig_regen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=scripts/figures/logs/r_fig_regen_%j.out
#SBATCH --error=scripts/figures/logs/r_fig_regen_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$(pwd)"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== Fig 2: Single-cell & Deconvolution (12 panels) ==="
Rscript scripts/figures/fig2_compact.R 2>&1 || echo "fig2_compact.R had errors (non-fatal)"

echo ""
echo "=== Fig 3 epigenomic panels ==="
Rscript scripts/figures/fig3_epigenomic_panels.R 2>&1 || echo "fig3_epigenomic_panels.R had errors (non-fatal)"

echo ""
echo "=== Supplementary ATAC R panels ==="
Rscript scripts/figures/figS_atac.R 2>&1 || echo "figS_atac.R had errors (non-fatal)"

echo ""
echo "=== Done: $(date) ==="
echo "Fig 2:"
ls -la figures/fig2/*.pdf 2>/dev/null || echo "  (none)"
echo "Fig 3 epigenomic:"
ls -la figures/fig3/panels/fig3_panel_f*.pdf figures/fig3/panels/fig3_panel_g*.pdf figures/fig3/panels/fig3_panel_h*.pdf 2>/dev/null || echo "  (none)"
echo "Supplementary ATAC:"
ls -la figures/supplementary/figS05_epigenomic_spatial/panel_09*.pdf figures/supplementary/figS05_epigenomic_spatial/panel_11*.pdf 2>/dev/null || echo "  (none)"
