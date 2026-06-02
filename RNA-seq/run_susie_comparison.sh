#!/bin/bash
# =============================================================================
# Launch post-hoc SuSiE comparison + cross-eQTL analysis
# Run AFTER SuSiE-COLOC and GTEx COLOC jobs complete
#
# Usage:
#   bash RNA-seq/run_susie_comparison.sh [SUSIE_JOB_IDS] [GTEX_JOB_ID]
#   e.g.: bash RNA-seq/run_susie_comparison.sh 14733346:14733347:14733348:14733349:14733350:14733351 14733328
# =============================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
LOG_DIR="${BASE}/RNA-seq/logs"

mkdir -p "${LOG_DIR}"

DEP_ARGS=""
if [[ "${1:-}" != "" ]] && [[ "${2:-}" != "" ]]; then
  DEP_ARGS="--dependency=afterok:${1},afterok:${2}"
  echo "Dependencies: SuSiE=${1}, GTEx=${2}"
elif [[ "${1:-}" != "" ]]; then
  DEP_ARGS="--dependency=afterok:${1}"
  echo "Dependencies: SuSiE=${1}"
fi

echo "=== Submitting comparison analyses ==="

# SuSiE vs ABF comparison
JOB1=$(sbatch ${DEP_ARGS} \
  --job-name="susie_abf_compare" \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=32G \
  --time=48:00:00 \
  --output="${LOG_DIR}/susie_abf_compare_%j.out" \
  --error="${LOG_DIR}/susie_abf_compare_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\"; micromamba activate rnaseq; cd ${BASE}; Rscript RNA-seq/35t_susie_abf_comparison.R'" \
  --parsable)
echo "  SuSiE vs ABF comparison: SLURM ${JOB1}"

# Cross-eQTL comparison
JOB2=$(sbatch ${DEP_ARGS} \
  --job-name="cross_eqtl_compare" \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=32G \
  --time=48:00:00 \
  --output="${LOG_DIR}/cross_eqtl_compare_%j.out" \
  --error="${LOG_DIR}/cross_eqtl_compare_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\"; micromamba activate rnaseq; cd ${BASE}; Rscript RNA-seq/35y_cross_eqtl_comparison.R'" \
  --parsable)
echo "  Cross-eQTL comparison: SLURM ${JOB2}"

echo "=== Comparison jobs submitted ==="
echo "After completion, run supplementary figures:"
echo "  sbatch scripts/figures/run_pub_figures.sh"
