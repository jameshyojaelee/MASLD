#!/bin/bash
#SBATCH --job-name=figS_dream_regen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/figS_dream_regen_%j.out
#SBATCH --error=scripts/figures/logs/figS_dream_regen_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

source ~/.bashrc
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

# Supplementary figures that load dream results and need ashr refresh
SCRIPTS=(
  "figS_sensitivity.R"
  "figS_deg_validation.R"
  "figS_combat_sensitivity.R"
  "figS_cohort_forest.R"
  "figS_sc_validation.R"
  "figS_mouse_overview.R"
)

for script_name in "${SCRIPTS[@]}"; do
  script="scripts/figures/${script_name}"
  if [[ -f "${script}" ]]; then
    run_figure "${script}"
  else
    echo "WARNING: ${script_name} not found, skipping"
  fi
done

echo "=========================================="
echo "SUPPLEMENTARY REGEN COMPLETE"
echo "=========================================="
if [[ ${#failed[@]} -gt 0 ]]; then
  echo "FAILED scripts (${#failed[@]}):"
  for f in "${failed[@]}"; do
    echo "  - ${f}"
  done
  exit 1
else
  echo "All supplementary figures regenerated successfully!"
  exit 0
fi
