#!/bin/bash
#SBATCH --job-name=compact_figs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/compact_figures_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/compact_figures_%j.err

set -euo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPT_DIR="${BASE}/scripts/figures"
FIG_DIR="${BASE}/figures"

mkdir -p "${FIG_DIR}/logs"

echo "=== Compact Publication Figures (4 Main) ==="
echo "Time: $(date)"
echo "Node: $(hostname)"
echo "BASE: ${BASE}"

# Activate environment
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

# Install ggrastr if missing (for rasterize_layer)
Rscript -e 'if (!requireNamespace("ggrastr", quietly=TRUE)) install.packages("ggrastr", repos="https://cloud.r-project.org")'

FAILED=0
SUCCESS=0

run_fig() {
  local name="$1"
  local script="$2"
  echo ""
  echo "--- ${name} ---"
  echo "Script: ${script}"
  echo "Start: $(date)"
  if Rscript "${script}"; then
    echo "  OK: ${name}"
    SUCCESS=$((SUCCESS + 1))
  else
    echo "  FAILED: ${name}"
    FAILED=$((FAILED + 1))
  fi
}

# --- Main Figures ---
run_fig "Fig 2: Single-cell & Deconvolution" "${SCRIPT_DIR}/fig2_compact.R"

# Sub-panels first (sourced by fig3_compact and fig4_compact)
run_fig "Fig 3 sub: Epigenomic"      "${SCRIPT_DIR}/fig3_epigenomic_panels.R"
run_fig "Fig 4 sub: Proteomics"      "${SCRIPT_DIR}/fig4_proteomics_panels.R"
run_fig "Fig 4 sub: Spatial"         "${SCRIPT_DIR}/fig4_spatial_panels.R"

# Composite figures (source sub-panel scripts)
run_fig "Fig 3: Causal & Epigenomic" "${SCRIPT_DIR}/fig3_compact.R"
run_fig "Fig 4: Pharma & Spatial"    "${SCRIPT_DIR}/fig4_compact.R"

echo ""
echo "=== Summary ==="
echo "Succeeded: ${SUCCESS}"
echo "Failed: ${FAILED}"
echo "Output: ${FIG_DIR}/"
echo "Done: $(date)"

if [ "${FAILED}" -gt 0 ]; then
  echo "WARNING: ${FAILED} figure(s) failed. Check logs above."
  exit 1
fi
