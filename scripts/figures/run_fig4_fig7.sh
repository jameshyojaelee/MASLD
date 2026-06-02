#!/bin/bash
#SBATCH --job-name=liver_fig4_fig7
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4_fig7_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4_fig7_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${BASE}/scripts/figures"

mkdir -p "${BASE}/figures" "${SCRIPT_DIR}/logs"
cd "${BASE}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

# Install ggrastr if missing
Rscript -e 'if (!requireNamespace("ggrastr", quietly=TRUE)) install.packages("ggrastr", repos="https://cloud.r-project.org", quiet=TRUE)' 2>/dev/null || true

export MASLD_PROJECT_ROOT="${BASE}"

echo "=== Figure 5: Causal Architecture ==="
echo "Start: $(date)"
Rscript "${SCRIPT_DIR}/fig5_causal_architecture.R"
echo "Fig 5 done: $(date)"
echo ""

echo "=== Figure 7: Multi-Evidence Atlas ==="
echo "Start: $(date)"
Rscript "${SCRIPT_DIR}/fig7_multi_evidence.R"
echo "Fig 7 done: $(date)"
echo ""

echo "=== Complete ==="
echo "End: $(date)"
echo ""
echo "Generated files:"
ls -lh "${BASE}/figures"/fig5*.pdf "${BASE}/figures"/fig7*.pdf 2>/dev/null || echo "  (no PDFs found)"
