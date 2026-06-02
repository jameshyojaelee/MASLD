#!/usr/bin/env bash
##############################################################################
# download_hra007511_aria2c.sh — Download HRA007511 Visium FASTQs via aria2c
#
# Usage: sbatch --export=BATCH_ID=1 download_hra007511_aria2c.sh
#   BATCH_ID: 1-8, each handles 5 runs (40 total runs / 8 batches)
#
# Features:
#   - HTTPS endpoint (faster than FTP)
#   - aria2c with 8 connections per file
#   - 2 files downloaded in parallel per batch
#   - Idempotent: skips completed files, resumes partial downloads
#   - Validates file size against server Content-Length
##############################################################################
#SBATCH --job-name=dl_hra_%x
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/data_retrieval/logs/dl_hra_batch%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/data_retrieval/logs/dl_hra_batch%x_%j.err

set -euo pipefail

###############################################################################
# Configuration
###############################################################################
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUTDIR="${PROJECT_ROOT}/data/HRA007511/fastq_visium"
HTTPS_BASE="https://download.cncb.ac.cn/gsa-human/HRA007511"
ARIA2C="/nfs/sw/easybuild/software/aria2/1.37.0-GCCcore-13.3.0/bin/aria2c"
LIBSSH2_DIR="/gpfs/commons/home/jameslee/scImmuneATLAS/.micromamba/pkgs/libssh2-1.11.1-h251f7ec_0/lib"
GCC_LIB="/nfs/sw/easybuild/software/GCCcore/13.3.0/lib64"

export LD_LIBRARY_PATH="${GCC_LIB}:${LIBSSH2_DIR}:${LD_LIBRARY_PATH:-}"

# aria2c settings
CONNECTIONS_PER_FILE=8   # -x and -s
MAX_CONCURRENT=2         # -j: parallel file downloads per batch
RETRY=10
RETRY_WAIT=30
TIMEOUT=120

###############################################################################
# All 40 runs (order matches spatial_datasets.yaml)
###############################################################################
ALL_RUNS=(
    # Batch 1: runs 1-5
    HRR1782834 HRR1782835 HRR1782836 HRR1782837 HRR1782838
    # Batch 2: runs 6-10
    HRR1782839 HRR1782840 HRR1782842 HRR1782843 HRR1782844
    # Batch 3: runs 11-15
    HRR1782845 HRR1782846 HRR1782847 HRR1782848 HRR1782849
    # Batch 4: runs 16-20
    HRR1782850 HRR1782851 HRR1782852 HRR1782853 HRR1782854
    # Batch 5: runs 21-25
    HRR2173353 HRR2173354 HRR2173355 HRR2173356 HRR2173357
    # Batch 6: runs 26-30
    HRR2173358 HRR2173359 HRR2173360 HRR2173361 HRR2173362
    # Batch 7: runs 31-35
    HRR2173363 HRR2173364 HRR2173365 HRR2173366 HRR2173367
    # Batch 8: runs 36-40
    HRR2173368 HRR2173369 HRR2173370 HRR2173371 HRR2173372
)

###############################################################################
# Determine batch slice
###############################################################################
BATCH_ID="${BATCH_ID:?ERROR: Set BATCH_ID=1..8 via --export}"
if (( BATCH_ID < 1 || BATCH_ID > 8 )); then
    echo "ERROR: BATCH_ID must be 1-8, got ${BATCH_ID}" >&2
    exit 1
fi

START_IDX=$(( (BATCH_ID - 1) * 5 ))
BATCH_RUNS=("${ALL_RUNS[@]:${START_IDX}:5}")

echo "=============================================="
echo "  HRA007511 Visium FASTQ Download"
echo "  Batch ${BATCH_ID}/8: ${BATCH_RUNS[*]}"
echo "  Node: $(hostname)"
echo "  Partition: ${SLURM_JOB_PARTITION:-unknown}"
echo "  Date: $(date)"
echo "  aria2c: ${ARIA2C}"
echo "=============================================="

# Verify aria2c works
if ! "${ARIA2C}" --version &>/dev/null; then
    echo "ERROR: aria2c not functional" >&2
    exit 1
fi

###############################################################################
# Build aria2c input file
###############################################################################
INPUT_FILE=$(mktemp /tmp/aria2_batch${BATCH_ID}_XXXXXX.txt)
trap "rm -f ${INPUT_FILE}" EXIT

SKIPPED=0
QUEUED=0

for RUN_ID in "${BATCH_RUNS[@]}"; do
    RUN_DIR="${OUTDIR}/${RUN_ID}"
    mkdir -p "${RUN_DIR}"

    for SUFFIX in "_f1.fastq.gz" "_r2.fastq.gz"; do
        FILENAME="${RUN_ID}${SUFFIX}"
        LOCAL_FILE="${RUN_DIR}/${FILENAME}"
        URL="${HTTPS_BASE}/${RUN_ID}/${FILENAME}"

        # Check if file already complete (compare to server Content-Length)
        if [ -f "${LOCAL_FILE}" ] && [ -s "${LOCAL_FILE}" ]; then
            LOCAL_SIZE=$(stat -c%s "${LOCAL_FILE}" 2>/dev/null || echo 0)
            REMOTE_SIZE=$(curl -skI "${URL}" 2>/dev/null | grep -i content-length | awk '{print $2}' | tr -d '\r' || echo 0)

            if [ "${LOCAL_SIZE}" = "${REMOTE_SIZE}" ] && [ "${REMOTE_SIZE}" != "0" ]; then
                echo "  SKIP: ${FILENAME} (${LOCAL_SIZE} bytes, matches server)"
                SKIPPED=$((SKIPPED + 1))
                continue
            else
                echo "  RESUME: ${FILENAME} (local=${LOCAL_SIZE}, server=${REMOTE_SIZE})"
            fi
        fi

        # Add to aria2c input file
        echo "${URL}"
        echo "  dir=${RUN_DIR}"
        echo "  out=${FILENAME}"
        echo ""
        QUEUED=$((QUEUED + 1))
    done
done >> "${INPUT_FILE}"

echo ""
echo "Skipped: ${SKIPPED} files (already complete)"
echo "Queued:  ${QUEUED} files for download"
echo ""

if [ "${QUEUED}" -eq 0 ]; then
    echo "All files already downloaded. Nothing to do."
    exit 0
fi

###############################################################################
# Run aria2c
###############################################################################
echo "Starting aria2c download (${MAX_CONCURRENT} parallel, ${CONNECTIONS_PER_FILE} connections/file)..."
echo "---"

"${ARIA2C}" \
    --input-file="${INPUT_FILE}" \
    --max-concurrent-downloads=${MAX_CONCURRENT} \
    --split=${CONNECTIONS_PER_FILE} \
    --max-connection-per-server=${CONNECTIONS_PER_FILE} \
    --min-split-size=10M \
    --max-tries=${RETRY} \
    --retry-wait=${RETRY_WAIT} \
    --timeout=${TIMEOUT} \
    --connect-timeout=60 \
    --check-certificate=false \
    --continue=true \
    --auto-file-renaming=false \
    --allow-overwrite=false \
    --summary-interval=60 \
    --console-log-level=notice \
    --file-allocation=none \
    --human-readable=true

EXIT_CODE=$?

echo ""
echo "=============================================="
echo "  Batch ${BATCH_ID} Complete"
echo "  Exit code: ${EXIT_CODE}"
echo "  Disk usage: $(du -sh "${OUTDIR}" 2>/dev/null | cut -f1)"
echo "  Date: $(date)"
echo "=============================================="

exit ${EXIT_CODE}
