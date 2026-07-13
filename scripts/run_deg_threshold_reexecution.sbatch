#!/bin/bash
#SBATCH --job-name=deg_threshold_reexec
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/deg_threshold_reexec_%j.log

# ===========================================================================
# DEG Threshold Migration Re-execution
# Regenerates all downstream outputs after threshold change:
#   OLD: lfsr < 0.05, |shrunk_logFC| > 0.2 → 7,588 DEGs
#   NEW: padj < 0.05, |logFC| > 0.3 → 5,484 DEGs
#
# Script 07 already re-run (5,484 DEGs confirmed). This script handles:
#   Step 1: Atlas rebuild (27a)
#   Step 2: Benchmark presets (27b)
#   Step 3: All figure scripts
#   Step 4: Supplementary sensitivity figures
#   Step 5: DEG threshold landscape figure
# ===========================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

mkdir -p logs

echo "=== DEG Threshold Re-execution started: $(date) ==="
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo ""

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

# ===================== Step 1: Atlas rebuild =====================
echo "====== STEP 1: Rebuild multi-evidence atlas ======"
run_script RNA-seq/27a_assemble_evidence_atlas.R

# ===================== Step 2: Benchmark presets =====================
echo "====== STEP 2: Benchmark presets ======"
run_script RNA-seq/27b_benchmark_presets.R

# ===================== Step 3: Core figure scripts =====================
echo "====== STEP 3: Figure scripts (DEG-dependent) ======"

# Main figures
run_script scripts/figures/fig5_translation.R
run_script scripts/figures/fig5_causal_architecture.R
# fig6_* scripts archived 2026-04-22 (6-fig → 5-fig restructure; content folded into fig5_convergence + figS_therapeutics)
# Legacy targets live in archive/fig6_ditched_2026-04-22/:
#   fig6_evidence_matrix_v2.R  fig6_upset_nonredundancy.R  fig6_convergence_sankey.R
#   fig6_therapeutic_windows.R fig6_pharmacotranscriptomics.R
run_script scripts/figures/fig5_convergence.R

# Supplementary figures
run_script scripts/figures/figS_sensitivity.R
run_script scripts/figures/figS_deg_threshold_landscape.R
run_script scripts/figures/figS_chen_influence.R
run_script scripts/figures/figS_depmap_essentiality.R
run_script scripts/figures/figS_ncrna.R
run_script scripts/figures/figS_cohort_forest.R
run_script scripts/figures/figS_sc_validation.R
run_script scripts/figures/fig_cohort_replication.R

# ===================== Step 4: Summary =====================
echo ""
echo "=========================================="
echo "=== DEG Threshold Re-execution Summary ==="
echo "=========================================="
echo "Succeeded (${#succeeded[@]}):"
for s in "${succeeded[@]}"; do echo "  ✓ $s"; done
echo ""
if [ ${#failed[@]} -gt 0 ]; then
  echo "FAILED (${#failed[@]}):"
  for f in "${failed[@]}"; do echo "  ✗ $f"; done
else
  echo "All scripts succeeded."
fi
echo ""
echo "=== Completed: $(date) ==="
