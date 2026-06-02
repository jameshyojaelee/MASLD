#!/bin/bash
# =============================================================================
# Launch SuSiE-COLOC array jobs for all 6 EUR GWAS x Broadaway eQTLs
# Each GWAS gets 22 array tasks (one per chromosome)
#
# Usage:
#   bash RNA-seq/run_susie_coloc_array.sh
# =============================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SBATCH_SCRIPT="${BASE}/RNA-seq/sbatch_susie_coloc_array.sh"
LOG_DIR="${BASE}/RNA-seq/logs"

mkdir -p "${LOG_DIR}"

GWAS_LIST=(UKBB_ALT UKBB_AST UKBB_GGT FINNGEN_NAFLD FINNGEN_NASH FINNGEN_HCC)

echo "=== Submitting SuSiE-COLOC array jobs (22 chr per GWAS) ==="
COMBINE_DEPS=""

for GWAS_NAME in "${GWAS_LIST[@]}"; do
  JOB_ID=$(sbatch \
    --job-name="susie_${GWAS_NAME}" \
    --output="${LOG_DIR}/susie_${GWAS_NAME}_chr%a_%A.out" \
    --error="${LOG_DIR}/susie_${GWAS_NAME}_chr%a_%A.err" \
    --export=ALL,GWAS_NAME="${GWAS_NAME}" \
    --parsable \
    "${SBATCH_SCRIPT}")
  echo "  ${GWAS_NAME}: array job ${JOB_ID} (22 tasks)"
  COMBINE_DEPS="${COMBINE_DEPS}:${JOB_ID}"
done

# Submit combiner job that runs after ALL array jobs complete
COMBINE_JOB=$(sbatch \
  --dependency="afterany${COMBINE_DEPS}" \
  --job-name="susie_combine" \
  --partition=cpu \
  --cpus-per-task=2 \
  --mem=16G \
  --time=48:00:00 \
  --output="${LOG_DIR}/susie_combine_%j.out" \
  --error="${LOG_DIR}/susie_combine_%j.err" \
  --wrap="cd ${BASE} && eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rnaseq && Rscript RNA-seq/35s_combine_susie_results.R" \
  --parsable)

echo ""
echo "  Combiner job: ${COMBINE_JOB} (runs after all array jobs)"
echo ""
echo "=== Total: ${#GWAS_LIST[@]} GWAS x 22 chr = $((${#GWAS_LIST[@]} * 22)) array tasks + 1 combiner ==="
