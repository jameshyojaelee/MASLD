#!/bin/bash
#SBATCH --job-name=fig3_supp7
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=scripts/figures/logs/fig3_supp7_%j.out
#SBATCH --error=scripts/figures/logs/fig3_supp7_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${BASE}/scripts/figures"
LOG_DIR="${SCRIPT_DIR}/logs"

mkdir -p "${BASE}/figures" "${LOG_DIR}"

cd "${BASE}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

echo "=========================================="
echo "[$(date '+%H:%M:%S')] Running fig3_deconvolution.R..."
echo "=========================================="
Rscript "${SCRIPT_DIR}/fig3_deconvolution.R" 2>&1
echo "[$(date '+%H:%M:%S')] fig3 done"
echo ""

echo "=========================================="
echo "[$(date '+%H:%M:%S')] Running supp_fig7_extended_celltype.R..."
echo "=========================================="
Rscript "${SCRIPT_DIR}/supp_fig7_extended_celltype.R" 2>&1
echo "[$(date '+%H:%M:%S')] supp_fig7 done"
echo ""

echo "=========================================="
echo "FIGURE GENERATION COMPLETE"
echo "=========================================="
ls -lh "${BASE}/figures/fig3_deconvolution.pdf" 2>/dev/null || echo "fig3 PDF not found"
ls -lh "${BASE}/figures/supp_fig7_extended_celltype.pdf" 2>/dev/null || echo "supp_fig7 PDF not found"
