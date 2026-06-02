#!/bin/bash
#SBATCH --job-name=ora_stg115_fix
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/ora_stg115_fix_%j.log

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
export MASLD_PROJECT_ROOT="${BASE}"

echo "=== Fix re-run: $(date) ==="

failed=()
succeeded=()

run_step() {
  local name="$1"
  shift
  echo "=========================================="
  echo "[$(date '+%H:%M:%S')] Running ${name}..."
  echo "=========================================="
  if "$@" 2>&1; then
    echo "[$(date '+%H:%M:%S')] ${name} SUCCEEDED"
    succeeded+=("${name}")
  else
    echo "[$(date '+%H:%M:%S')] ${name} FAILED"
    failed+=("${name}")
  fi
  echo ""
}

run_step "02c_ora_concordance" Rscript Analysis/Cross_Species_Concordance/scripts/02c_ora_concordance.R
run_step "115_transition_programs" Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/115_transition_programs.R

echo "=========================================="
echo "Succeeded (${#succeeded[@]}):"
for s in "${succeeded[@]}"; do echo "  ✓ $s"; done
if [ ${#failed[@]} -gt 0 ]; then
  echo "FAILED (${#failed[@]}):"
  for f in "${failed[@]}"; do echo "  ✗ $f"; done
fi
echo "=== Completed: $(date) ==="
