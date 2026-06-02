#!/bin/bash
#SBATCH --job-name=atlas_rebuild_web
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/atlas_rebuild_web_%j.log

# Rebuild multi-evidence atlas + web app data after DEG threshold change

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
export MASLD_PROJECT_ROOT="${BASE}"

echo "=== Atlas rebuild + web preprocessing: $(date) ==="

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

# Step 1: Rebuild multi-evidence atlas
run_step "27a_assemble_evidence_atlas" Rscript RNA-seq/27a_assemble_evidence_atlas.R

# Step 2: Benchmark presets (validates atlas integrity)
run_step "27b_benchmark_presets" Rscript RNA-seq/27b_benchmark_presets.R

# Step 3: Web app data preprocessing
run_step "preprocess_atlas_data" python masld-atlas-v2/scripts/preprocess_atlas_data.py

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
