#!/usr/bin/env bash
#SBATCH --job-name=dl_hra_visium
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/download_hra007511_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/download_hra007511_%j.err
##############################################################################
# 00a_download_hra007511_visium.sh — Download HRA007511 Visium FASTQs from NGDC
#
# 35 Visium experiments, 40 runs (paired-end: _f1 + _r2)
# Total: 80 FASTQ files, estimated ~200-300 GB
#
# Uses wget with retry logic. Can be re-run safely (skips existing files).
##############################################################################
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUTDIR="data/HRA007511/fastq_visium"
mkdir -p "${OUTDIR}"

FTP_BASE="ftp://download.big.ac.cn/gsa-human/HRA007511"

# Batch 1: 20 runs from experiments 11-25
BATCH1_RUNS=(
    HRR1782834 HRR1782835 HRR1782836 HRR1782837 HRR1782838
    HRR1782839 HRR1782840 HRR1782842 HRR1782843 HRR1782844
    HRR1782845 HRR1782846 HRR1782847 HRR1782848 HRR1782849
    HRR1782850 HRR1782851 HRR1782852 HRR1782853 HRR1782854
)

# Batch 2: 20 runs from experiments 27-29, 45-61
BATCH2_RUNS=(
    HRR2173353 HRR2173354 HRR2173355 HRR2173356 HRR2173357
    HRR2173358 HRR2173359 HRR2173360 HRR2173361 HRR2173362
    HRR2173363 HRR2173364 HRR2173365 HRR2173366 HRR2173367
    HRR2173368 HRR2173369 HRR2173370 HRR2173371 HRR2173372
)

ALL_RUNS=("${BATCH1_RUNS[@]}" "${BATCH2_RUNS[@]}")

echo "=============================================="
echo "  HRA007511 Visium FASTQ Download"
echo "  ${#ALL_RUNS[@]} runs × 2 files = $((${#ALL_RUNS[@]} * 2)) FASTQs"
echo "  Date: $(date)"
echo "=============================================="

TOTAL=${#ALL_RUNS[@]}
COMPLETED=0
FAILED=0

for RUN_ID in "${ALL_RUNS[@]}"; do
    RUN_DIR="${OUTDIR}/${RUN_ID}"
    mkdir -p "${RUN_DIR}"

    F1="${RUN_DIR}/${RUN_ID}_f1.fastq.gz"
    R2="${RUN_DIR}/${RUN_ID}_r2.fastq.gz"

    COMPLETED=$((COMPLETED + 1))
    echo ""
    echo "[${COMPLETED}/${TOTAL}] ${RUN_ID}"

    # Download R1 (_f1)
    if [ -f "${F1}" ] && [ -s "${F1}" ]; then
        echo "  R1: exists ($(du -h "${F1}" | cut -f1))"
    else
        echo "  R1: downloading..."
        wget -q --tries=5 --retry-connrefused --timeout=60 \
            -O "${F1}" \
            "${FTP_BASE}/${RUN_ID}/${RUN_ID}_f1.fastq.gz" \
            && echo "  R1: done ($(du -h "${F1}" | cut -f1))" \
            || { echo "  R1: FAILED"; FAILED=$((FAILED + 1)); rm -f "${F1}"; }
    fi

    # Download R2 (_r2)
    if [ -f "${R2}" ] && [ -s "${R2}" ]; then
        echo "  R2: exists ($(du -h "${R2}" | cut -f1))"
    else
        echo "  R2: downloading..."
        wget -q --tries=5 --retry-connrefused --timeout=60 \
            -O "${R2}" \
            "${FTP_BASE}/${RUN_ID}/${RUN_ID}_r2.fastq.gz" \
            && echo "  R2: done ($(du -h "${R2}" | cut -f1))" \
            || { echo "  R2: FAILED"; FAILED=$((FAILED + 1)); rm -f "${R2}"; }
    fi
done

echo ""
echo "=============================================="
echo "  Download Complete"
echo "  Total runs: ${TOTAL}"
echo "  Failed files: ${FAILED}"
echo "  Output: ${OUTDIR}/"
echo "  Disk usage: $(du -sh "${OUTDIR}" | cut -f1)"
echo "  Date: $(date)"
echo "=============================================="

if [ "${FAILED}" -gt 0 ]; then
    echo "WARNING: ${FAILED} files failed. Re-run this script to retry."
    exit 1
fi
