#!/bin/bash
#SBATCH --job-name=coloc_figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/coloc_figures_%j.out
#SBATCH --error=scripts/figures/logs/coloc_figures_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${BASE}/scripts/figures"

cd "${BASE}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript -e 'if (!requireNamespace("ggrastr", quietly=TRUE)) install.packages("ggrastr", repos="https://cloud.r-project.org", quiet=TRUE)' 2>/dev/null || true

export MASLD_PROJECT_ROOT="${BASE}"

failed=()

run_figure() {
  local script="$1"
  local name=$(basename "$script" .R)
  echo "=========================================="
  echo "[$(date '+%H:%M:%S')] Running ${name}..."
  echo "=========================================="
  if Rscript "${script}"; then
    echo "[$(date '+%H:%M:%S')] ${name} SUCCEEDED"
  else
    echo "[$(date '+%H:%M:%S')] ${name} FAILED"
    failed+=("${name}")
  fi
  echo ""
}

# Run the 3 COLOC-related figure scripts
run_figure "${SCRIPT_DIR}/fig3_compact.R"
run_figure "${SCRIPT_DIR}/figS_causal_extended.R"
run_figure "${SCRIPT_DIR}/figS_coloc_sensitivity.R"
run_figure "${SCRIPT_DIR}/figS_sensitivity.R"

echo "=========================================="
echo "COLOC FIGURE GENERATION COMPLETE"
echo "=========================================="

if [[ ${#failed[@]} -gt 0 ]]; then
  echo "FAILED scripts (${#failed[@]}):"
  for f in "${failed[@]}"; do
    echo "  - ${f}"
  done
  exit 1
else
  echo "All COLOC figures generated successfully!"
  exit 0
fi
