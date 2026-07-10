#!/bin/bash
#SBATCH --job-name=figregen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/rerender_3fixes_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/rerender_3fixes_%j.err
set -o pipefail
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export MASLD_PROJECT_ROOT="$BASE"; cd "$BASE"
eval "$(micromamba shell hook --shell bash)"; micromamba activate rnaseq
PASS=(); FAIL=()
run(){ local k="$1"; shift; echo "=== [$k] $* ==="; local t0=$SECONDS
  if [ "$k" = R ]; then timeout 1200 Rscript "scripts/figures/$1"; else timeout 1200 python "scripts/figures/$@"; fi
  local rc=$?; local dt=$((SECONDS-t0))
  if [ $rc -eq 0 ]; then PASS+=("$* (${dt}s)"); echo "[PASS] $* (${dt}s)"; else FAIL+=("$* rc=$rc"); echo "[FAIL] $* rc=$rc (${dt}s)"; fi; }
run PY published_deg_comparison.py
run R  figS_batch_correction.R
echo ""; echo "#### 3-FIX SUMMARY: ${#PASS[@]} PASS / ${#FAIL[@]} FAIL ####"
echo "PASS:"; printf '  %s\n' "${PASS[@]:-none}"; echo "FAIL:"; printf '  %s\n' "${FAIL[@]:-none}"
