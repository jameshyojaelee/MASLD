#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=phase10_followup
#SBATCH --output=scripts/figures/logs/phase10_followup_%j.out
#SBATCH --error=scripts/figures/logs/phase10_followup_%j.err

# Phase 10 Team C followup: assemble fig3 + re-render the 4-way EUR LD-panel
# comparison figures. These need to run after the array (assemble_fig3) or
# after Team A's path-rewritten archive scripts are in place.

set -eo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${PROJECT_ROOT}"

source ~/.bashrc
micromamba activate rnaseq

run_step () {
  local script="$1"
  echo "============================================================"
  echo "Running: ${script}"
  echo "Started: $(date)"
  echo "============================================================"
  if Rscript "${script}"; then
    echo "[PASS] ${script}"
  else
    echo "[FAIL] ${script} (continuing with remaining steps)"
  fi
}

run_step "scripts/figures/assemble_fig3.R"
run_step "scripts/figures/presentation_4way_eur.R"
run_step "scripts/figures/presentation_ld_panel_comparison.R"
run_step "scripts/figures/figS_polyfun_fm_concordance.R"

echo "============================================================"
echo "Followup complete: $(date)"
echo "============================================================"
