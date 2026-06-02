#!/bin/bash
#SBATCH --job-name=regen_fig5_fig6_fix
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/regen_fig5_fig6_fix_%j.out
#SBATCH --error=scripts/figures/logs/regen_fig5_fig6_fix_%j.err

# Re-run 3 figure scripts that failed with duplicate factor / missing column errors

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${BASE}/scripts/figures"
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

run_figure "${SCRIPT_DIR}/fig5_causal_architecture.R"
run_figure "${SCRIPT_DIR}/fig5_translation.R"
# fig6_therapeutic_windows.R — archived 2026-04-22 (6-fig → 5-fig restructure; content folded into fig5_convergence)

echo "=========================================="
if [[ ${#failed[@]} -gt 0 ]]; then
  echo "FAILED (${#failed[@]}):"
  for f in "${failed[@]}"; do echo "  - ${f}"; done
  exit 1
else
  echo "All 3 figures regenerated successfully!"
fi
