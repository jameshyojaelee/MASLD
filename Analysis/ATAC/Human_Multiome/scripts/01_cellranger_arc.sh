#!/bin/bash
#SBATCH --job-name=cellranger_arc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --array=1-18
#SBATCH --output=logs/cellranger_arc_%A_%a.out
#SBATCH --error=logs/cellranger_arc_%A_%a.err

# =============================================================================
# CellRanger ARC count for GSE244832 10x Multiome (snRNA + snATAC)
# Module 2, Step 2a: Process paired RNA + ATAC FASTQs per donor
#
# Each array task processes one donor from sample_pairing.csv:
#   - Creates a libraries.csv with RNA and ATAC FASTQ paths
#   - Handles non-standard FASTQ naming (symlinks to 10x convention)
#   - Runs cellranger-arc count
#
# Prerequisites:
#   1. Fill in rna_srrs in sample_pairing.csv (see comment header for donor SRR lists)
#   2. Ensure FASTQs are downloaded to data/GSE244832/fastq/{SRR_ID}/
#
# Usage:
#   cd Analysis/ATAC/Human_Multiome
#   mkdir -p logs
#   sbatch scripts/01_cellranger_arc.sh
# =============================================================================

set -euo pipefail

# --- Configuration -----------------------------------------------------------
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FASTQ_BASE="${PROJECT_ROOT}/data/GSE244832/fastq"
SCRIPT_DIR="${PROJECT_ROOT}/Analysis/ATAC/Human_Multiome/scripts"
OUTPUT_BASE="${PROJECT_ROOT}/Analysis/ATAC/Human_Multiome/cellranger_arc"
ARC_REF="/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A"
PAIRING_CSV="${SCRIPT_DIR}/sample_pairing.csv"

LOCALCORES=16
LOCALMEM=120

# --- Load module --------------------------------------------------------------
module load cellranger-arc/2.0.2

# --- Parse sample pairing CSV ------------------------------------------------
# Skip comment lines (starting with #) and header line, then get the Nth data row
TASK_ID=${SLURM_ARRAY_TASK_ID}
LINE=$(grep -v '^#' "${PAIRING_CSV}" | tail -n +2 | sed -n "${TASK_ID}p")

if [[ -z "${LINE}" ]]; then
    echo "ERROR: No data found for array task ${TASK_ID} in ${PAIRING_CSV}"
    exit 1
fi

DONOR_ID=$(echo "${LINE}" | cut -d',' -f1)
ATAC_SRR=$(echo "${LINE}" | cut -d',' -f2)
RNA_SRRS=$(echo "${LINE}" | cut -d',' -f3)

echo "============================================================"
echo "CellRanger ARC count"
echo "  Task ID:    ${TASK_ID}"
echo "  Donor:      ${DONOR_ID}"
echo "  ATAC SRR:   ${ATAC_SRR}"
echo "  RNA SRRs:   ${RNA_SRRS}"
echo "  Started:    $(date)"
echo "============================================================"

# --- Validate inputs ----------------------------------------------------------
if [[ "${RNA_SRRS}" == "PLACEHOLDER" ]]; then
    echo "ERROR: RNA SRRs not filled in for ${DONOR_ID}. Edit sample_pairing.csv first."
    exit 1
fi

# Function to ensure FASTQ directory has 10x-compatible naming
# CellRanger ARC requires: {sample}_S{N}_L{NNN}_R{1,2,3}_001.fastq.gz
# SRA downloads sometimes produce: {SRR}_1.fastq.gz, {SRR}_2.fastq.gz
ensure_10x_naming() {
    local fastq_dir="$1"
    local srr_id="$2"

    if [[ ! -d "${fastq_dir}" ]]; then
        echo "ERROR: FASTQ directory does not exist: ${fastq_dir}"
        return 1
    fi

    # Check if already 10x-compatible
    if ls "${fastq_dir}"/${srr_id}_S1_L001_R1_001.fastq.gz &>/dev/null; then
        return 0
    fi

    # Check for non-standard naming and create symlinks
    if ls "${fastq_dir}"/${srr_id}_1.fastq.gz &>/dev/null; then
        echo "  Creating 10x-compatible symlinks in ${fastq_dir}"
        for f in "${fastq_dir}"/${srr_id}_*.fastq.gz; do
            base=$(basename "$f")
            # Extract read number: SRR_1.fastq.gz -> 1
            read_num=$(echo "$base" | sed "s/${srr_id}_\([0-9]*\)\.fastq\.gz/\1/")
            target="${fastq_dir}/${srr_id}_S1_L001_R${read_num}_001.fastq.gz"
            if [[ ! -e "${target}" ]]; then
                ln -s "$(basename "$f")" "${target}"
                echo "    Symlinked: $(basename "$f") -> $(basename "${target}")"
            fi
        done
        return 0
    fi

    echo "ERROR: No recognized FASTQ files in ${fastq_dir} for ${srr_id}"
    ls "${fastq_dir}"/
    return 1
}

# Validate and fix ATAC FASTQ directory
ATAC_FASTQ_DIR="${FASTQ_BASE}/${ATAC_SRR}"
echo ""
echo "Validating ATAC FASTQs..."
ensure_10x_naming "${ATAC_FASTQ_DIR}" "${ATAC_SRR}"

# Validate and fix RNA FASTQ directories (semicolon-separated SRR list)
echo ""
echo "Validating RNA FASTQs..."
IFS=';' read -ra RNA_SRR_ARRAY <<< "${RNA_SRRS}"
for rna_srr in "${RNA_SRR_ARRAY[@]}"; do
    rna_fastq_dir="${FASTQ_BASE}/${rna_srr}"
    ensure_10x_naming "${rna_fastq_dir}" "${rna_srr}"
done

# --- Build libraries.csv -----------------------------------------------------
WORK_DIR="${OUTPUT_BASE}/${DONOR_ID}"
mkdir -p "${WORK_DIR}"

LIBRARIES_CSV="${WORK_DIR}/libraries.csv"
echo "fastqs,sample,library_type" > "${LIBRARIES_CSV}"

# Add ATAC FASTQ directory
echo "${ATAC_FASTQ_DIR},${ATAC_SRR},Chromatin Accessibility" >> "${LIBRARIES_CSV}"

# Add RNA FASTQ directories (one row per SRR)
for rna_srr in "${RNA_SRR_ARRAY[@]}"; do
    rna_fastq_dir="${FASTQ_BASE}/${rna_srr}"
    echo "${rna_fastq_dir},${rna_srr},Gene Expression" >> "${LIBRARIES_CSV}"
done

echo ""
echo "Generated libraries.csv:"
cat "${LIBRARIES_CSV}"
echo ""

# --- Run CellRanger ARC count ------------------------------------------------
# CellRanger ARC creates output in the current working directory under --id
cd "${OUTPUT_BASE}"

echo "Running cellranger-arc count..."
echo "  Output dir: ${OUTPUT_BASE}/${DONOR_ID}"
echo ""

cellranger-arc count \
    --id="${DONOR_ID}" \
    --reference="${ARC_REF}" \
    --libraries="${LIBRARIES_CSV}" \
    --localcores=${LOCALCORES} \
    --localmem=${LOCALMEM}

echo ""
echo "============================================================"
echo "CellRanger ARC count COMPLETE for ${DONOR_ID}"
echo "  Finished:   $(date)"
echo "  Output:     ${OUTPUT_BASE}/${DONOR_ID}/outs/"
echo "============================================================"
