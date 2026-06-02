#!/usr/bin/env bash
##############################################################################
# submit_hra007511_all_batches.sh — Submit all 8 download batches
#
# Spreads jobs across io (2 nodes) and cpu (many nodes) partitions:
#   - Batches 1-4: io partition (2 per node)
#   - Batches 5-8: cpu partition (separate nodes)
#
# Each batch downloads 5 runs (10 FASTQs) with aria2c multi-connection.
# All 8 jobs run simultaneously = 8 × 2 parallel files = 16 files at once.
#
# Usage: bash submit_hra007511_all_batches.sh
##############################################################################
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOWNLOAD_SCRIPT="${SCRIPT_DIR}/download_hra007511_aria2c.sh"
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

echo "=============================================="
echo "  Submitting 8 parallel HRA007511 download jobs"
echo "  Date: $(date)"
echo "=============================================="

declare -a JOB_IDS=()

# Batches 1-4 → io partition (optimized for I/O, 2 nodes available)
for BATCH in 1 2 3 4; do
    JOB_ID=$(sbatch --parsable \
        --job-name="dl_hra_b${BATCH}" \
        --partition=io \
        --cpus-per-task=4 \
        --mem=8G \
        --time=48:00:00 \
        --output="${LOG_DIR}/dl_hra_batch${BATCH}_%j.out" \
        --error="${LOG_DIR}/dl_hra_batch${BATCH}_%j.err" \
        --export="BATCH_ID=${BATCH}" \
        "${DOWNLOAD_SCRIPT}")
    JOB_IDS+=("${JOB_ID}")
    echo "  Batch ${BATCH}/8 → io partition  → Job ${JOB_ID}"
done

# Batches 5-8 → cpu partition (many nodes available)
for BATCH in 5 6 7 8; do
    JOB_ID=$(sbatch --parsable \
        --job-name="dl_hra_b${BATCH}" \
        --partition=cpu \
        --cpus-per-task=4 \
        --mem=8G \
        --time=48:00:00 \
        --output="${LOG_DIR}/dl_hra_batch${BATCH}_%j.out" \
        --error="${LOG_DIR}/dl_hra_batch${BATCH}_%j.err" \
        --export="BATCH_ID=${BATCH}" \
        "${DOWNLOAD_SCRIPT}")
    JOB_IDS+=("${JOB_ID}")
    echo "  Batch ${BATCH}/8 → cpu partition → Job ${JOB_ID}"
done

echo ""
echo "=============================================="
echo "  All 8 jobs submitted!"
echo "  Job IDs: ${JOB_IDS[*]}"
echo ""
echo "  Monitor: squeue -u \$USER -n dl_hra_b1,dl_hra_b2,dl_hra_b3,dl_hra_b4,dl_hra_b5,dl_hra_b6,dl_hra_b7,dl_hra_b8"
echo "  Logs:    ${LOG_DIR}/dl_hra_batch*"
echo ""
echo "  Expected: ~1.2 TB across 40 runs (80 FASTQs)"
echo "  Output:   data/HRA007511/fastq_visium/"
echo "=============================================="

# Save job IDs for monitoring
echo "${JOB_IDS[*]}" > "${LOG_DIR}/hra007511_job_ids.txt"
echo "Job IDs saved to ${LOG_DIR}/hra007511_job_ids.txt"
