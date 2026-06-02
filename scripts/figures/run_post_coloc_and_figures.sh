#!/bin/bash
#SBATCH --job-name=post_coloc_figs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/post_coloc_figs_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/post_coloc_figs_%j.err

# Post-COLOC pipeline: rebuild atlas, cross-ancestry, benchmark, then all figures
# Submit with: sbatch --dependency=afterok:14525579:14525580:14525581:14525693:14525694:14525695 scripts/figures/run_post_coloc_and_figures.sh

set -euo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPT_DIR="${BASE}/scripts/figures"
FIG_DIR="${BASE}/figures"

mkdir -p "${FIG_DIR}/logs"

echo "=== Post-COLOC Pipeline + All Publication Figures ==="
echo "Time: $(date)"
echo "Node: $(hostname)"
echo "BASE: ${BASE}"

# Activate environment
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

# Install ggrastr if missing
Rscript -e 'if (!requireNamespace("ggrastr", quietly=TRUE)) install.packages("ggrastr", repos="https://cloud.r-project.org")'

FAILED=0
SUCCESS=0

run_step() {
  local name="$1"
  local cmd="$2"
  echo ""
  echo "--- ${name} ---"
  echo "Start: $(date)"
  if eval "${cmd}"; then
    echo "  OK: ${name} ($(date))"
    SUCCESS=$((SUCCESS + 1))
  else
    echo "  FAILED: ${name} ($(date))"
    FAILED=$((FAILED + 1))
  fi
}

# =========================================================================
# Phase 1: Post-COLOC analyses
# =========================================================================

echo ""
echo "=== Phase 1: Post-COLOC Analyses ==="

# Verify COLOC outputs exist
echo "Checking COLOC outputs..."
# ARCHIVED 2026-04-08: FinnGen removed (finngen_nafld finngen_nash finngen_hcc)
for d in bbj_alt bbj_ast bbj_ggt; do
  f="${BASE}/RNA-seq/results/causal_inference/${d}/coloc_results.csv"
  if [ -f "$f" ]; then
    echo "  OK: ${d}/coloc_results.csv"
  else
    echo "  MISSING: ${d}/coloc_results.csv"
  fi
done

# Step 3a: Rebuild multi-evidence atlas
run_step "Atlas rebuild (27a)" "cd ${BASE}/RNA-seq && Rscript 27a_assemble_evidence_atlas.R"

# Step 3b: Cross-ancestry replication
run_step "Cross-ancestry (48)" "cd ${BASE}/RNA-seq && Rscript 48_cross_ancestry_replication.R"

# Step 3c: Benchmark presets
run_step "Benchmark (27b)" "cd ${BASE}/RNA-seq && Rscript 27b_benchmark_presets.R"

# =========================================================================
# Phase 2: All publication figures
# =========================================================================

echo ""
echo "=== Phase 2: Publication Figures ==="

run_step "Fig 1: Atlas"               "Rscript ${SCRIPT_DIR}/fig1_compact.R"
run_step "Fig 2: Single-cell & Deconvolution" "Rscript ${SCRIPT_DIR}/fig2_compact.R"

# Sub-panels (sourced by fig3/fig4 compact)
run_step "Fig 3 sub: Epigenomic"      "Rscript ${SCRIPT_DIR}/fig3_epigenomic_panels.R"
run_step "Fig 4 sub: Proteomics"      "Rscript ${SCRIPT_DIR}/fig4_proteomics_panels.R"
run_step "Fig 4 sub: Spatial"         "Rscript ${SCRIPT_DIR}/fig4_spatial_panels.R"

# Composite figures
run_step "Fig 3: Causal & Epigenomic" "Rscript ${SCRIPT_DIR}/fig3_compact.R"
run_step "Fig 4: Pharma & Spatial"    "Rscript ${SCRIPT_DIR}/fig4_compact.R"

# =========================================================================
# Summary
# =========================================================================

echo ""
echo "=== Summary ==="
echo "Succeeded: ${SUCCESS}"
echo "Failed: ${FAILED}"
echo "Output: ${FIG_DIR}/"
echo ""

# List output PDFs
echo "Generated PDFs:"
ls -lh "${FIG_DIR}"/*.pdf 2>/dev/null || echo "  (none)"

echo ""
echo "Done: $(date)"

if [ "${FAILED}" -gt 0 ]; then
  echo "WARNING: ${FAILED} step(s) failed. Check logs above."
  exit 1
fi
