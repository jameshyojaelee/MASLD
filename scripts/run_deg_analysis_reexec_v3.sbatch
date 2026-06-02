#!/bin/bash
#SBATCH --job-name=deg_analysis_v3
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/deg_analysis_v3_%j.log

# Re-run 2 scripts that failed with 'tier' column error after DEG schema change

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
export MASLD_PROJECT_ROOT="${BASE}"

echo "=== Analysis re-execution v3 started: $(date) ==="

failed=()
succeeded=()

run_script() {
  local script="$1"
  local name=$(basename "$script" .R)
  echo "=========================================="
  echo "[$(date '+%H:%M:%S')] Running ${name}..."
  echo "=========================================="
  if Rscript "${script}" 2>&1; then
    echo "[$(date '+%H:%M:%S')] ${name} SUCCEEDED"
    succeeded+=("${name}")
  else
    echo "[$(date '+%H:%M:%S')] ${name} FAILED (continuing)"
    failed+=("${name}")
  fi
  echo ""
}

run_script RNA-seq/20_drug_repurposing_v3.R
run_script RNA-seq/29_positive_control_validation.R

echo ""
echo "=========================================="
echo "=== Summary ==="
echo "=========================================="
echo "Succeeded (${#succeeded[@]}):"
for s in "${succeeded[@]}"; do echo "  ✓ $s"; done
if [ ${#failed[@]} -gt 0 ]; then
  echo "FAILED (${#failed[@]}):"
  for f in "${failed[@]}"; do echo "  ✗ $f"; done
fi
echo "=== Completed: $(date) ==="
