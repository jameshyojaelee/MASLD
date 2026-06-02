#!/bin/bash
#SBATCH --job-name=cellranger_atac
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --array=1-18
#SBATCH --output=logs/cellranger_atac_%A_%a.out
#SBATCH --error=logs/cellranger_atac_%A_%a.err

# =============================================================================
# CellRanger-ATAC count for GSE244832 snATAC-seq (18 donors)
#
# IMPORTANT: GSE244832 is NOT a true 10x Multiome. The snRNA-seq used 10x
# Chromium 3' v3.1 and snATAC-seq used combinatorial barcoding on separate
# nuclei. Therefore we use CellRanger-ATAC (not CellRanger ARC).
#
# Cell type labels will be transferred computationally from the separately
# processed snRNA-seq data in downstream SnapATAC2 processing (script 02).
#
# Usage:
#   cd Analysis/ATAC/Human_Multiome
#   mkdir -p logs
#   sbatch scripts/01_cellranger_atac.sh
# =============================================================================

set -euo pipefail

# --- Configuration ---
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FASTQ_BASE="${PROJECT_ROOT}/data/GSE244832/fastq"
OUTPUT_BASE="${PROJECT_ROOT}/Analysis/ATAC/Human_Multiome/cellranger_atac"
ATAC_REF="/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-atac-GRCh38-2024-A"
DONOR_PAIRING="${PROJECT_ROOT}/data/GSE244832/metadata/donor_pairing.csv"

LOCALCORES=16
LOCALMEM=120

# --- Load module ---
module load cellranger-atac/2.1.0

# --- Parse donor pairing CSV ---
TASK_ID=${SLURM_ARRAY_TASK_ID}
# Skip header, get Nth data row
LINE=$(tail -n +2 "${DONOR_PAIRING}" | sed -n "${TASK_ID}p")

if [[ -z "${LINE}" ]]; then
    echo "ERROR: No data found for array task ${TASK_ID}"
    exit 1
fi

DONOR_ID=$(echo "${LINE}" | cut -d',' -f1)
ATAC_SAMPLE=$(echo "${LINE}" | cut -d',' -f2)
ATAC_SRR=$(echo "${LINE}" | cut -d',' -f3)
CONDITION=$(echo "${LINE}" | cut -d',' -f8)

echo "============================================================"
echo "CellRanger-ATAC count"
echo "  Task ID:    ${TASK_ID}"
echo "  Donor:      ${DONOR_ID} (${ATAC_SAMPLE})"
echo "  ATAC SRR:   ${ATAC_SRR}"
echo "  Condition:  ${CONDITION}"
echo "  Started:    $(date)"
echo "  Node:       $(hostname)"
echo "============================================================"

# --- Validate FASTQ directory ---
ATAC_FASTQ_DIR="${FASTQ_BASE}/${ATAC_SRR}"
if [[ ! -d "${ATAC_FASTQ_DIR}" ]]; then
    echo "ERROR: FASTQ directory not found: ${ATAC_FASTQ_DIR}"
    exit 1
fi

# Ensure 10x-compatible naming (SRA downloads may use non-standard names)
if ls "${ATAC_FASTQ_DIR}"/${ATAC_SRR}_S1_L001_R1_001.fastq.gz &>/dev/null; then
    echo "FASTQs: 10x naming OK"
elif ls "${ATAC_FASTQ_DIR}"/${ATAC_SRR}_1.fastq.gz &>/dev/null; then
    echo "Creating 10x-compatible symlinks..."
    for f in "${ATAC_FASTQ_DIR}"/${ATAC_SRR}_*.fastq.gz; do
        base=$(basename "$f")
        read_num=$(echo "$base" | sed "s/${ATAC_SRR}_\([0-9]*\)\.fastq\.gz/\1/")
        target="${ATAC_FASTQ_DIR}/${ATAC_SRR}_S1_L001_R${read_num}_001.fastq.gz"
        if [[ ! -e "${target}" ]]; then
            ln -s "$(basename "$f")" "${target}"
            echo "  Symlinked: $(basename "$f") -> $(basename "${target}")"
        fi
    done
else
    echo "ERROR: No recognized FASTQ files in ${ATAC_FASTQ_DIR}"
    ls "${ATAC_FASTQ_DIR}"/
    exit 1
fi

echo ""
echo "FASTQ files:"
ls -lh "${ATAC_FASTQ_DIR}"/*.fastq.gz 2>/dev/null | head -10
echo ""

# --- Validate reference ---
if [[ ! -d "${ATAC_REF}" ]]; then
    echo "ERROR: CellRanger-ATAC reference not found: ${ATAC_REF}"
    exit 1
fi

# --- Run CellRanger-ATAC count ---
mkdir -p "${OUTPUT_BASE}"
cd "${OUTPUT_BASE}"

# Remove incomplete previous run if exists
if [[ -d "${DONOR_ID}" ]] && [[ ! -f "${DONOR_ID}/outs/fragments.tsv.gz" ]]; then
    echo "Removing incomplete previous run for ${DONOR_ID}..."
    rm -rf "${DONOR_ID}"
fi

# Skip if already completed
if [[ -f "${DONOR_ID}/outs/fragments.tsv.gz" ]]; then
    echo "Output already exists for ${DONOR_ID}, skipping."
    exit 0
fi

echo "Running cellranger-atac count..."
cellranger-atac count \
    --id="${DONOR_ID}" \
    --reference="${ATAC_REF}" \
    --fastqs="${ATAC_FASTQ_DIR}" \
    --sample="${ATAC_SRR}" \
    --localcores=${LOCALCORES} \
    --localmem=${LOCALMEM}

echo ""
echo "============================================================"
echo "CellRanger-ATAC COMPLETE for ${DONOR_ID} (${ATAC_SAMPLE})"
echo "  Finished:   $(date)"
echo "  Fragments:  ${OUTPUT_BASE}/${DONOR_ID}/outs/fragments.tsv.gz"
echo "  Peaks:      ${OUTPUT_BASE}/${DONOR_ID}/outs/peaks.bed"
echo "============================================================"
