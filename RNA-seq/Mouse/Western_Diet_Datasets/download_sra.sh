#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=logs/sra_download_%A_%a.out
#SBATCH --error=logs/sra_download_%A_%a.err

# =============================================================================
# Download raw FASTQs from SRA for Western Diet mouse datasets
# =============================================================================
# Two-step process: prefetch (resolve + cache) -> fasterq-dump (extract FASTQ)
# This is more reliable than fasterq-dump alone, which can fail to resolve
# older accessions.
#
# Usage:
#   cd RNA-seq/Mouse/Western_Diet_Datasets
#   sbatch --array=1-17%10 download_sra.sh GSE220575
#   sbatch --array=1-30%10 download_sra.sh GSE246088
#   sbatch --array=1-39%10 download_sra.sh GSE305484
#
# Prerequisites:
#   - <DATASET_ID>/metadata/srr_list.txt with one SRR accession per line
#
# Output:
#   <DATASET_ID>/fastq/<SRR>_1.fastq.gz  (R1, paired-end)
#   <DATASET_ID>/fastq/<SRR>_2.fastq.gz  (R2, paired-end)
#   <DATASET_ID>/fastq/<SRR>.fastq.gz    (single-end, if applicable)
# =============================================================================

set -euo pipefail

module load SRA-Toolkit/3.2.0-gompi-2023b
module load pigz/2.7-GCCcore-12.2.0

# ---- Configuration ----
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets"
DATASET_ID="${1:?ERROR: Supply DATASET_ID as first argument (e.g., GSE220575)}"
SRR_LIST="${BASE}/${DATASET_ID}/metadata/srr_list.txt"
OUTDIR="${BASE}/${DATASET_ID}/fastq"
# Use scratch for temp and prefetch cache (avoid filling home directory)
TMPDIR_LOCAL="/scratch/jameslee/sra_tmp/${DATASET_ID}/${SLURM_ARRAY_TASK_ID}"
PREFETCH_DIR="/scratch/jameslee/sra_prefetch/${DATASET_ID}"

# ---- Validate ----
if [[ ! -f "${SRR_LIST}" ]]; then
    echo "ERROR: SRR list not found: ${SRR_LIST}"
    exit 1
fi

TOTAL=$(wc -l < "${SRR_LIST}")
if [[ ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL} ]]; then
    echo "SKIP: Task ${SLURM_ARRAY_TASK_ID} > total SRRs (${TOTAL})"
    exit 0
fi

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SRR_LIST}")
if [[ -z "${SRR}" ]]; then
    echo "ERROR: Empty SRR at line ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

mkdir -p "${OUTDIR}" "${TMPDIR_LOCAL}" "${PREFETCH_DIR}"

echo "============================================="
echo "  SRA FASTQ Download (prefetch + fasterq-dump)"
echo "============================================="
echo "  Dataset:   ${DATASET_ID}"
echo "  SRR:       ${SRR} (task ${SLURM_ARRAY_TASK_ID}/${TOTAL})"
echo "  Output:    ${OUTDIR}"
echo "  Temp:      ${TMPDIR_LOCAL}"
echo "  Prefetch:  ${PREFETCH_DIR}"
echo "  Started:   $(date)"
echo "============================================="

# ---- Check if already downloaded ----
if [[ -f "${OUTDIR}/${SRR}_1.fastq.gz" && -f "${OUTDIR}/${SRR}_2.fastq.gz" ]]; then
    R1_SIZE=$(stat -c%s "${OUTDIR}/${SRR}_1.fastq.gz" 2>/dev/null || echo 0)
    R2_SIZE=$(stat -c%s "${OUTDIR}/${SRR}_2.fastq.gz" 2>/dev/null || echo 0)
    if [[ ${R1_SIZE} -gt 1000000 && ${R2_SIZE} -gt 1000000 ]]; then
        echo "SKIP: FASTQs already exist and look complete (R1=${R1_SIZE}, R2=${R2_SIZE} bytes)"
        exit 0
    fi
fi
if [[ -f "${OUTDIR}/${SRR}.fastq.gz" ]]; then
    SE_SIZE=$(stat -c%s "${OUTDIR}/${SRR}.fastq.gz" 2>/dev/null || echo 0)
    if [[ ${SE_SIZE} -gt 1000000 ]]; then
        echo "SKIP: Single-end FASTQ already exists (${SE_SIZE} bytes)"
        exit 0
    fi
fi

START_TIME=$SECONDS

# ---- Step 1: Prefetch SRA file ----
echo "[$(date)] Step 1: prefetch ${SRR}..."
prefetch "${SRR}" \
    --max-size 50G \
    --output-directory "${PREFETCH_DIR}" \
    --force ALL \
    2>&1

ELAPSED_PREFETCH=$(( SECONDS - START_TIME ))
echo "[$(date)] Prefetch completed in ${ELAPSED_PREFETCH}s"

# Find the prefetched .sra file
SRA_FILE=$(find "${PREFETCH_DIR}/${SRR}" -name "*.sra" -o -name "${SRR}" 2>/dev/null | head -1)
if [[ -z "${SRA_FILE}" || ! -f "${SRA_FILE}" ]]; then
    # Fallback: check the default cache location
    SRA_FILE="${PREFETCH_DIR}/${SRR}/${SRR}.sra"
    if [[ ! -f "${SRA_FILE}" ]]; then
        echo "ERROR: Prefetched file not found in ${PREFETCH_DIR}/${SRR}/"
        ls -la "${PREFETCH_DIR}/${SRR}/" 2>/dev/null
        exit 1
    fi
fi
echo "  SRA file: ${SRA_FILE} ($(du -h "${SRA_FILE}" | cut -f1))"

# ---- Step 2: fasterq-dump from cached file ----
echo "[$(date)] Step 2: fasterq-dump (extract FASTQ)..."
START_FQD=$SECONDS

fasterq-dump \
    --split-files \
    --threads ${SLURM_CPUS_PER_TASK} \
    --outdir "${TMPDIR_LOCAL}" \
    --temp "${TMPDIR_LOCAL}" \
    "${SRA_FILE}"

ELAPSED_FQD=$(( SECONDS - START_FQD ))
echo "[$(date)] fasterq-dump completed in ${ELAPSED_FQD}s"

# ---- Step 3: Compress with pigz ----
echo "[$(date)] Step 3: Compressing with pigz..."
START_COMPRESS=$SECONDS

if [[ -f "${TMPDIR_LOCAL}/${SRR}_1.fastq" && -f "${TMPDIR_LOCAL}/${SRR}_2.fastq" ]]; then
    pigz -p ${SLURM_CPUS_PER_TASK} "${TMPDIR_LOCAL}/${SRR}_1.fastq"
    pigz -p ${SLURM_CPUS_PER_TASK} "${TMPDIR_LOCAL}/${SRR}_2.fastq"
    mv "${TMPDIR_LOCAL}/${SRR}_1.fastq.gz" "${OUTDIR}/"
    mv "${TMPDIR_LOCAL}/${SRR}_2.fastq.gz" "${OUTDIR}/"
    echo "  Paired-end: ${SRR}_1.fastq.gz, ${SRR}_2.fastq.gz"
elif [[ -f "${TMPDIR_LOCAL}/${SRR}.fastq" ]]; then
    pigz -p ${SLURM_CPUS_PER_TASK} "${TMPDIR_LOCAL}/${SRR}.fastq"
    mv "${TMPDIR_LOCAL}/${SRR}.fastq.gz" "${OUTDIR}/"
    echo "  Single-end: ${SRR}.fastq.gz"
else
    echo "ERROR: No FASTQ files found after fasterq-dump"
    ls -la "${TMPDIR_LOCAL}/" | grep -i "${SRR}"
    exit 1
fi

ELAPSED_COMPRESS=$(( SECONDS - START_COMPRESS ))

# ---- Step 4: Clean up prefetch cache for this accession ----
rm -rf "${PREFETCH_DIR}/${SRR}" 2>/dev/null
rm -rf "${TMPDIR_LOCAL}" 2>/dev/null
echo "[$(date)] Cleaned up temp files"

# ---- Report ----
echo ""
echo "============================================="
echo "  Download Summary"
echo "============================================="
echo "  SRR:         ${SRR}"
echo "  Prefetch:    ${ELAPSED_PREFETCH}s"
echo "  fasterq-dump: ${ELAPSED_FQD}s"
echo "  Compress:    ${ELAPSED_COMPRESS}s"
echo "  Total:       $(( SECONDS - START_TIME ))s"
ls -lh "${OUTDIR}/${SRR}"*fastq.gz 2>/dev/null
echo "  Completed:   $(date)"
echo "============================================="
