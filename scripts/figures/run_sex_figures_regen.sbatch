#!/bin/bash
#SBATCH --job-name=sex_fig_regen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/sex_fig_regen_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/sex_fig_regen_%j.err

set -euo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPT_DIR="${BASE}/scripts/figures"

mkdir -p "${BASE}/figures/logs"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

# Install ggrastr if missing
Rscript -e 'if (!requireNamespace("ggrastr", quietly=TRUE)) install.packages("ggrastr", repos="https://cloud.r-project.org")' 2>/dev/null || true

echo "=== Regenerating sex-related figures ==="
echo "Date: $(date)"
echo "Node: $(hostname)"

echo ""
echo "--- fig2_disease_progression.R (Panels i,j: sex scatter + class counts) ---"
Rscript "${SCRIPT_DIR}/fig2_disease_progression.R"

echo ""
echo "--- figS_sex_divergent.R (Panels a-c: sex divergent sensitivity) ---"
Rscript "${SCRIPT_DIR}/figS_sex_divergent.R"

echo ""
echo "=== Done: $(date) ==="
echo "Outputs:"
ls -lh "${BASE}/figures/"*disease_progression* 2>/dev/null || true
ls -lh "${BASE}/figures/supplementary/figS03_deconvolution/figS_sex_divergent.pdf" 2>/dev/null || true
ls -lh "${BASE}/figures/main/fig2_f2_switch/fig2i_"* 2>/dev/null || true
ls -lh "${BASE}/figures/main/fig2_f2_switch/fig2j_"* 2>/dev/null || true
