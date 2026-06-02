#!/bin/bash
#SBATCH --job-name=fig345_regen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/fig345_regen_%j.out
#SBATCH --error=scripts/figures/logs/fig345_regen_%j.err

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
  echo "[$(date)] Running ${name}..."
  echo "=========================================="
  if Rscript "${script}"; then
    echo "[$(date)] ${name} SUCCEEDED"
  else
    echo "[$(date)] ${name} FAILED (exit $?)"
    failed+=("${name}")
  fi
  echo ""
}

# Figure 3: Disease Progression Atlas
run_figure scripts/figures/fig3_progression.R

# Figure 4: Pharma + Proteomics + Spatial (sources sub-scripts)
run_figure scripts/figures/fig4_compact.R

# Figure 5: Translational Convergence & Validation
run_figure scripts/figures/fig5_translation.R

echo "=========================================="
echo "[$(date)] Summary"
echo "=========================================="
echo "Outputs:"
ls -lh figures/main/fig3_deconfounding/fig3_progression.pdf 2>/dev/null || echo "  fig3: MISSING"
ls -lh figures/main/fig4_causal_architecture/fig4_compact.pdf 2>/dev/null || echo "  fig4: MISSING"
ls -lh figures/main/fig5_sex_dimorphism/fig5_translation.pdf 2>/dev/null || echo "  fig5: MISSING"

if [ ${#failed[@]} -gt 0 ]; then
  echo ""
  echo "FAILED scripts: ${failed[*]}"
  exit 1
else
  echo ""
  echo "All 3 figures regenerated successfully."
fi
