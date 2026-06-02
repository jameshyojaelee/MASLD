#!/bin/bash
# =============================================================================
# Launch SuSiE-COLOC for EUR GWAS x Broadaway eQTLs
#
# Usage:
#   bash RNA-seq/run_susie_coloc.sh [LD_PREP_JOB_ID]
#
# If LD_PREP_JOB_ID is provided, adds SLURM dependency (--dependency=afterok:ID)
# Note: FinnGen studies archived 2026-04-08 (probable duplicates).
#       Whitfield/2023_36653562_* studies archived 2026-04-09 (provenance unverified).
# =============================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SBATCH_SCRIPT="${BASE}/RNA-seq/sbatch_susie_coloc.sh"
LOG_DIR="${BASE}/RNA-seq/logs"

mkdir -p "${LOG_DIR}"

LD_DEP=""
if [[ "${1:-}" != "" ]]; then
  LD_DEP="--dependency=afterok:${1}"
  echo "LD prep dependency: job ${1}"
fi

# FinnGen studies removed 2026-04-08 (archived as probable duplicates)
# To add more GWAS, append to this list (must match gwas_registry.tsv study IDs)
GWAS_LIST=(UKBB_ALT UKBB_AST UKBB_GGT)

echo "=== Submitting SuSiE-COLOC jobs ==="
for GWAS_NAME in "${GWAS_LIST[@]}"; do
  JOB_ID=$(sbatch ${LD_DEP} \
    --job-name="susie_${GWAS_NAME}" \
    --output="${LOG_DIR}/susie_${GWAS_NAME}_%j.out" \
    --error="${LOG_DIR}/susie_${GWAS_NAME}_%j.err" \
    --export=ALL,GWAS_NAME="${GWAS_NAME}" \
    --parsable \
    "${SBATCH_SCRIPT}")
  echo "  ${GWAS_NAME}: SLURM ${JOB_ID}"
done

echo "=== All SuSiE-COLOC jobs submitted ==="
