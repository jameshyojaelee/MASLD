#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --output=logs/sra_batch_%A_%a.out
#SBATCH --error=logs/sra_batch_%A_%a.err

# =============================================================================
# Batch SRA download: processes ALL samples for one dataset sequentially
# =============================================================================
# Since interactive QOS limits to 4 concurrent jobs, each job processes
# all samples in one dataset sequentially. Use 3 jobs (one per dataset).
#
# Usage (submit all 3 datasets as array tasks):
#   sbatch --array=1-3 download_sra_batch.sh
#
# Array mapping:
#   1 = GSE220575 (17 samples)
#   2 = GSE246088 (30 samples)
#   3 = GSE305484 (39 samples)
# =============================================================================

set -euo pipefail

module load SRA-Toolkit/3.2.0-gompi-2023b
module load pigz/2.7-GCCcore-12.2.0

# ---- Configuration ----
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets"

# Map array task to dataset
case ${SLURM_ARRAY_TASK_ID} in
    1) DATASET_ID="GSE220575" ;;
    2) DATASET_ID="GSE246088" ;;
    3) DATASET_ID="GSE305484" ;;
    *) echo "ERROR: Invalid task ID ${SLURM_ARRAY_TASK_ID}"; exit 1 ;;
esac

SRR_LIST="${BASE}/${DATASET_ID}/metadata/srr_list.txt"
OUTDIR="${BASE}/${DATASET_ID}/fastq"
PREFETCH_DIR="/scratch/jameslee/sra_prefetch/${DATASET_ID}"
TMPDIR_BASE="/scratch/jameslee/sra_tmp/${DATASET_ID}"

# ---- Validate ----
if [[ ! -f "${SRR_LIST}" ]]; then
    echo "ERROR: SRR list not found: ${SRR_LIST}"
    exit 1
fi

TOTAL=$(wc -l < "${SRR_LIST}")
mkdir -p "${OUTDIR}" "${PREFETCH_DIR}"

echo "============================================="
echo "  Batch SRA Download"
echo "============================================="
echo "  Dataset:   ${DATASET_ID}"
echo "  Samples:   ${TOTAL}"
echo "  Output:    ${OUTDIR}"
echo "  Started:   $(date)"
echo "============================================="
echo ""

COMPLETED=0
SKIPPED=0
FAILED=0
TOTAL_START=$SECONDS

while IFS= read -r SRR; do
    [[ -z "${SRR}" ]] && continue
    IDX=$((COMPLETED + SKIPPED + FAILED + 1))
    TMPDIR="${TMPDIR_BASE}/${SRR}"
    mkdir -p "${TMPDIR}"

    echo "--- [${IDX}/${TOTAL}] ${SRR} ---"

    # ---- Check if already downloaded ----
    if [[ -f "${OUTDIR}/${SRR}_1.fastq.gz" && -f "${OUTDIR}/${SRR}_2.fastq.gz" ]]; then
        R1_SIZE=$(stat -c%s "${OUTDIR}/${SRR}_1.fastq.gz" 2>/dev/null || echo 0)
        R2_SIZE=$(stat -c%s "${OUTDIR}/${SRR}_2.fastq.gz" 2>/dev/null || echo 0)
        if [[ ${R1_SIZE} -gt 1000000 && ${R2_SIZE} -gt 1000000 ]]; then
            echo "  SKIP: FASTQs exist (R1=${R1_SIZE}, R2=${R2_SIZE})"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi
    fi
    if [[ -f "${OUTDIR}/${SRR}.fastq.gz" ]]; then
        SE_SIZE=$(stat -c%s "${OUTDIR}/${SRR}.fastq.gz" 2>/dev/null || echo 0)
        if [[ ${SE_SIZE} -gt 1000000 ]]; then
            echo "  SKIP: SE FASTQ exists (${SE_SIZE})"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi
    fi

    SAMPLE_START=$SECONDS

    # ---- Step 1: Prefetch ----
    echo "  [$(date +%H:%M:%S)] Prefetch..."
    if ! prefetch "${SRR}" --max-size 50G --output-directory "${PREFETCH_DIR}" --force ALL 2>&1; then
        echo "  ERROR: prefetch failed for ${SRR}"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Find prefetched file
    SRA_FILE=$(find "${PREFETCH_DIR}/${SRR}" -name "*.sra" -o -name "${SRR}" 2>/dev/null | head -1)
    if [[ -z "${SRA_FILE}" || ! -f "${SRA_FILE}" ]]; then
        SRA_FILE="${PREFETCH_DIR}/${SRR}/${SRR}.sra"
    fi
    if [[ ! -f "${SRA_FILE}" ]]; then
        echo "  ERROR: Prefetched file not found for ${SRR}"
        FAILED=$((FAILED + 1))
        continue
    fi
    echo "  SRA: $(du -h "${SRA_FILE}" | cut -f1)"

    # ---- Step 2: fasterq-dump ----
    echo "  [$(date +%H:%M:%S)] fasterq-dump..."
    if ! fasterq-dump --split-files --threads ${SLURM_CPUS_PER_TASK} --outdir "${TMPDIR}" --temp "${TMPDIR}" "${SRA_FILE}" 2>&1; then
        echo "  ERROR: fasterq-dump failed for ${SRR}"
        FAILED=$((FAILED + 1))
        rm -rf "${PREFETCH_DIR}/${SRR}" "${TMPDIR}" 2>/dev/null
        continue
    fi

    # ---- Step 3: Compress ----
    echo "  [$(date +%H:%M:%S)] Compressing..."
    if [[ -f "${TMPDIR}/${SRR}_1.fastq" && -f "${TMPDIR}/${SRR}_2.fastq" ]]; then
        pigz -p ${SLURM_CPUS_PER_TASK} "${TMPDIR}/${SRR}_1.fastq"
        pigz -p ${SLURM_CPUS_PER_TASK} "${TMPDIR}/${SRR}_2.fastq"
        mv "${TMPDIR}/${SRR}_1.fastq.gz" "${OUTDIR}/"
        mv "${TMPDIR}/${SRR}_2.fastq.gz" "${OUTDIR}/"
        R1_FINAL=$(du -h "${OUTDIR}/${SRR}_1.fastq.gz" | cut -f1)
        R2_FINAL=$(du -h "${OUTDIR}/${SRR}_2.fastq.gz" | cut -f1)
        echo "  PE: R1=${R1_FINAL}, R2=${R2_FINAL}"
    elif [[ -f "${TMPDIR}/${SRR}.fastq" ]]; then
        pigz -p ${SLURM_CPUS_PER_TASK} "${TMPDIR}/${SRR}.fastq"
        mv "${TMPDIR}/${SRR}.fastq.gz" "${OUTDIR}/"
        SE_FINAL=$(du -h "${OUTDIR}/${SRR}.fastq.gz" | cut -f1)
        echo "  SE: ${SE_FINAL}"
    else
        echo "  ERROR: No FASTQ output for ${SRR}"
        FAILED=$((FAILED + 1))
        rm -rf "${PREFETCH_DIR}/${SRR}" "${TMPDIR}" 2>/dev/null
        continue
    fi

    # Cleanup
    rm -rf "${PREFETCH_DIR}/${SRR}" "${TMPDIR}" 2>/dev/null
    ELAPSED=$((SECONDS - SAMPLE_START))
    COMPLETED=$((COMPLETED + 1))
    echo "  Done in ${ELAPSED}s (${COMPLETED} done, ${SKIPPED} skipped, ${FAILED} failed)"
    echo ""

done < "${SRR_LIST}"

TOTAL_ELAPSED=$(( SECONDS - TOTAL_START ))
echo "============================================="
echo "  ${DATASET_ID} Download Summary"
echo "============================================="
echo "  Total samples:  ${TOTAL}"
echo "  Completed:      ${COMPLETED}"
echo "  Skipped:        ${SKIPPED}"
echo "  Failed:         ${FAILED}"
echo "  Wall time:      ${TOTAL_ELAPSED}s ($(( TOTAL_ELAPSED / 3600 ))h $(( (TOTAL_ELAPSED % 3600) / 60 ))m)"
echo "  Finished:       $(date)"
echo "============================================="

# Final disk usage
echo "  Disk usage:"
du -sh "${OUTDIR}" 2>/dev/null
echo ""

if [[ ${FAILED} -gt 0 ]]; then
    echo "WARNING: ${FAILED} samples failed. Check logs above for details."
    exit 1
fi
