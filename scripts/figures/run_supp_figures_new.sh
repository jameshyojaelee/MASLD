#!/bin/bash
#SBATCH --job-name=supp_figs_new
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=scripts/figures/logs/supp_figs_new_%j.out
#SBATCH --error=scripts/figures/logs/supp_figs_new_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${BASE}/scripts/figures"
FIG_OUT="${BASE}/figures"
LOG_DIR="${SCRIPT_DIR}/logs"

mkdir -p "${FIG_OUT}" "${LOG_DIR}"

cd "${BASE}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

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

# New supplementary figures
NEW_FIGURES=(
  "figS_cohort_forest.R"
  "figS_mouse_overview.R"
  "figS_evidence_convergence.R"
  "figS_sc_validation.R"
)

for script_name in "${NEW_FIGURES[@]}"; do
  script="${SCRIPT_DIR}/${script_name}"
  if [[ -f "${script}" ]]; then
    run_figure "${script}"
  else
    echo "WARNING: ${script_name} not found"
  fi
done

echo "=========================================="
echo "SUPPLEMENTARY FIGURE GENERATION COMPLETE"
echo "=========================================="
echo ""
echo "Generated files:"
ls -lh "${FIG_OUT}"/figS_*.pdf 2>/dev/null || echo "  (no PDFs found)"
echo ""

if [[ ${#failed[@]} -gt 0 ]]; then
  echo "FAILED scripts (${#failed[@]}):"
  for f in "${failed[@]}"; do
    echo "  - ${f}"
  done
  exit 1
else
  echo "All supplementary figures generated successfully!"
  exit 0
fi
