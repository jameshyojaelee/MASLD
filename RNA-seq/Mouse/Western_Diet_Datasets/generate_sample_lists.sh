#!/bin/bash
# =============================================================================
# Generate sample_list.txt for Kallisto quant from downloaded FASTQs
# =============================================================================
# Run AFTER download_sra.sh completes. Generates the tab-separated sample_list.txt
# needed by run_kallisto_quant.sh.
#
# Usage:
#   bash generate_sample_lists.sh GSE220575
#   bash generate_sample_lists.sh GSE246088
#   bash generate_sample_lists.sh GSE305484
#   bash generate_sample_lists.sh all    # generate for all 3
# =============================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets"

generate_list() {
    local DATASET_ID="$1"
    local SRR_LIST="${BASE}/${DATASET_ID}/metadata/srr_list.txt"
    local FASTQDIR="${BASE}/${DATASET_ID}/fastq"
    local OUTFILE="${BASE}/${DATASET_ID}/metadata/sample_list.txt"

    echo "=== Generating sample_list.txt for ${DATASET_ID} ==="

    if [[ ! -f "${SRR_LIST}" ]]; then
        echo "ERROR: SRR list not found: ${SRR_LIST}"
        return 1
    fi

    > "${OUTFILE}"  # truncate

    local TOTAL=0
    local FOUND_PE=0
    local FOUND_SE=0
    local MISSING=0

    while IFS= read -r SRR; do
        [[ -z "${SRR}" ]] && continue
        TOTAL=$((TOTAL + 1))

        R1="${FASTQDIR}/${SRR}_1.fastq.gz"
        R2="${FASTQDIR}/${SRR}_2.fastq.gz"
        SE="${FASTQDIR}/${SRR}.fastq.gz"

        if [[ -f "${R1}" && -f "${R2}" ]]; then
            printf '%s\t%s\t%s\n' "${SRR}" "${R1}" "${R2}" >> "${OUTFILE}"
            FOUND_PE=$((FOUND_PE + 1))
        elif [[ -f "${SE}" ]]; then
            printf '%s\t%s\n' "${SRR}" "${SE}" >> "${OUTFILE}"
            FOUND_SE=$((FOUND_SE + 1))
        else
            echo "  WARNING: No FASTQ for ${SRR}"
            MISSING=$((MISSING + 1))
        fi
    done < "${SRR_LIST}"

    echo "  Total SRRs:    ${TOTAL}"
    echo "  Paired-end:    ${FOUND_PE}"
    echo "  Single-end:    ${FOUND_SE}"
    echo "  Missing:       ${MISSING}"
    echo "  Output:        ${OUTFILE}"
    echo ""
}

# Parse argument
DATASET="${1:?Usage: bash generate_sample_lists.sh <DATASET_ID|all>}"

if [[ "${DATASET}" == "all" ]]; then
    for ds in GSE220575 GSE246088 GSE305484; do
        generate_list "${ds}"
    done
else
    generate_list "${DATASET}"
fi
