#!/bin/bash
# =============================================================================
# Launch GTEx v8 Liver × UKBB ALT replication COLOC
# =============================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SBATCH_SCRIPT="${BASE}/RNA-seq/sbatch_gtex_coloc.sh"
LOG_DIR="${BASE}/RNA-seq/logs"

mkdir -p "${LOG_DIR}"

echo "=== Submitting GTEx replication COLOC ==="
JOB_ID=$(sbatch \
  --output="${LOG_DIR}/gtex_coloc_ukbb_%j.out" \
  --error="${LOG_DIR}/gtex_coloc_ukbb_%j.err" \
  --parsable \
  "${SBATCH_SCRIPT}")
echo "  GTEx × UKBB ALT: SLURM ${JOB_ID}"
echo "=== Done ==="
