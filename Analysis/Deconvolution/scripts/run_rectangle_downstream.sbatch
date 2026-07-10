#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/rectangle_downstream_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/rectangle_downstream_%j.err
#
# Rectangle downstream: runs AFTER the deconvolution array (18756501) writes all
# 5 cohorts' *_rectangle_proportions.tsv. Submit with:
#   sbatch --dependency=afterok:<ARRAY_JOBID> run_rectangle_downstream.sbatch
#
# Steps (each non-fatal to the next; the heavy DE runs first):
#   1. 25b_deconv_attribution_rectangle.R  -> RNA-seq/results/causal_inference/rectangle/   (limma-voom-qw, composition-adjusted DE; the Fig3-parallel payload)
#   2. 29_compare_three_methods.R          -> Analysis/Deconvolution/rectangle_comparison/  (Rectangle vs InstaPrism/BayesPrism vs MuSiC concordance)
#   3. 30_attribution_method_robustness.R  -> Analysis/Deconvolution/rectangle_comparison/  (Jaccard of hepatocyte-intrinsic sets across methods; needs step 1)
# All figures obey the project rules; NO existing file is modified (new dirs only).
# ---------------------------------------------------------------------------

set -uo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="${ROOT}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "[$(date)] Rectangle downstream on $(hostname), CPUs=${SLURM_CPUS_PER_TASK:-NA}"

rc=0
run_step () {  # $1 = label, $2 = Rscript path
  echo ""
  echo "========== [$(date)] STEP: $1 =========="
  if Rscript "$2"; then
    echo "[$(date)] STEP OK: $1"
  else
    echo "[$(date)] STEP FAILED (continuing): $1" >&2
    rc=1
  fi
}

# 1. Heavy composition-adjusted DE (Rectangle proportions swapped in) — runs first.
run_step "25b attribution (Rectangle)" \
  "${ROOT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/25b_deconv_attribution_rectangle.R"

# 2. 3-way proportion concordance + Unknown-content summary.
run_step "29 compare three methods" \
  "${ROOT}/Analysis/Deconvolution/scripts/29_compare_three_methods.R"

# 3. Attribution-class robustness (needs step 1 output).
run_step "30 attribution robustness" \
  "${ROOT}/Analysis/Deconvolution/scripts/30_attribution_method_robustness.R"

echo ""
echo "[$(date)] Rectangle downstream complete (aggregate rc=${rc})."
exit ${rc}
