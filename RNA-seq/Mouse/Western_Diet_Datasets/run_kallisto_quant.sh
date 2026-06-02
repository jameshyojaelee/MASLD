#!/bin/bash
#SBATCH --job-name=kallisto_quant
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=logs/kallisto_quant_%A_%a.out
#SBATCH --error=logs/kallisto_quant_%A_%a.err

# =============================================================================
# Kallisto quantification pipeline for Western Diet mouse datasets
# =============================================================================
# Usage:
#   # Paired-end (default):
#   sbatch --array=1-$(wc -l < sample_list.txt) run_kallisto_quant.sh <DATASET_ID>
#
#   # Single-end (supply fragment length + SD):
#   sbatch --array=1-$(wc -l < sample_list.txt) run_kallisto_quant.sh <DATASET_ID> --single -l 200 -s 30
#
# Input:
#   <DATASET_ID>/metadata/sample_list.txt — tab-separated file:
#     Paired-end: SRR_accession \t R1_path \t R2_path
#     Single-end: SRR_accession \t fastq_path
#
# Output:
#   <DATASET_ID>/quant/kallisto/<SRR>/abundance.h5
#   <DATASET_ID>/quant/kallisto/<SRR>/abundance.tsv
#   <DATASET_ID>/quant/kallisto/<SRR>/run_info.json
#
# Requires: module load kallisto/0.51.1
# =============================================================================

set -euo pipefail

module load kallisto/0.51.1

# ---- Configuration ----
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INDEX="${BASE}/data/reference/kallisto/gencode.vM38.kallisto.idx"
GTF="${BASE}/data/ncrna_conservation/gencode.vM38.annotation.gtf.gz"
WDDIR="${BASE}/RNA-seq/Mouse/Western_Diet_Datasets"

# ---- Parse arguments ----
DATASET_ID="${1:?ERROR: Must supply DATASET_ID as first argument (e.g., GSE220575)}"
shift  # Remaining args are extra kallisto flags (e.g., --single -l 200 -s 30)
EXTRA_FLAGS=("$@")

SAMPLE_LIST="${WDDIR}/${DATASET_ID}/metadata/sample_list.txt"
OUTBASE="${WDDIR}/${DATASET_ID}/quant/kallisto"

# ---- Validate ----
if [[ ! -f "${INDEX}" ]]; then
    echo "ERROR: Kallisto index not found: ${INDEX}"
    echo "Run build_kallisto_index.sh first."
    exit 1
fi

if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: Sample list not found: ${SAMPLE_LIST}"
    echo "Expected format (tab-separated):"
    echo "  Paired-end: SRR_accession <TAB> R1_path <TAB> R2_path"
    echo "  Single-end: SRR_accession <TAB> fastq_path"
    exit 1
fi

# ---- Extract sample info for this array task ----
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")
if [[ -z "${LINE}" ]]; then
    echo "ERROR: No entry for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

NCOLS=$(echo "${LINE}" | awk -F'\t' '{print NF}')
SRR=$(echo "${LINE}" | cut -f1)
OUTDIR="${OUTBASE}/${SRR}"
mkdir -p "${OUTDIR}"

echo "============================================="
echo "  Kallisto Quantification"
echo "============================================="
echo "  Dataset:  ${DATASET_ID}"
echo "  Sample:   ${SRR} (task ${SLURM_ARRAY_TASK_ID})"
echo "  Index:    ${INDEX}"
echo "  Threads:  ${SLURM_CPUS_PER_TASK}"
echo "  Output:   ${OUTDIR}"

# ---- Run Kallisto ----
START_TIME=$SECONDS

if [[ ${NCOLS} -ge 3 ]]; then
    # Paired-end
    R1=$(echo "${LINE}" | cut -f2)
    R2=$(echo "${LINE}" | cut -f3)
    echo "  Mode:     paired-end"
    echo "  R1:       ${R1}"
    echo "  R2:       ${R2}"
    echo "============================================="

    # Validate FASTQ existence
    for f in "${R1}" "${R2}"; do
        if [[ ! -f "${f}" ]]; then
            echo "ERROR: FASTQ not found: ${f}"
            exit 1
        fi
    done

    kallisto quant \
        -i "${INDEX}" \
        -o "${OUTDIR}" \
        -t "${SLURM_CPUS_PER_TASK}" \
        "${EXTRA_FLAGS[@]+"${EXTRA_FLAGS[@]}"}" \
        "${R1}" "${R2}"

elif [[ ${NCOLS} -eq 2 ]]; then
    # Single-end — requires --single -l <frag_len> -s <frag_sd> in EXTRA_FLAGS
    FQ=$(echo "${LINE}" | cut -f2)
    echo "  Mode:     single-end"
    echo "  FASTQ:    ${FQ}"
    echo "  Extra:    ${EXTRA_FLAGS[*]+"${EXTRA_FLAGS[*]}"}"
    echo "============================================="

    if [[ ! -f "${FQ}" ]]; then
        echo "ERROR: FASTQ not found: ${FQ}"
        exit 1
    fi

    # Ensure --single is passed for SE data
    HAS_SINGLE=0
    for flag in "${EXTRA_FLAGS[@]+"${EXTRA_FLAGS[@]}"}"; do
        [[ "${flag}" == "--single" ]] && HAS_SINGLE=1
    done
    if [[ ${HAS_SINGLE} -eq 0 ]]; then
        echo "WARNING: Single-end FASTQ detected but --single not in extra flags."
        echo "  Adding --single -l 200 -s 30 (defaults). Override with explicit args."
        EXTRA_FLAGS+=(--single -l 200 -s 30)
    fi

    kallisto quant \
        -i "${INDEX}" \
        -o "${OUTDIR}" \
        -t "${SLURM_CPUS_PER_TASK}" \
        "${EXTRA_FLAGS[@]}" \
        "${FQ}"

else
    echo "ERROR: Unexpected column count (${NCOLS}) in sample list line ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

ELAPSED=$(( SECONDS - START_TIME ))

# ---- Verify output ----
if [[ -f "${OUTDIR}/abundance.tsv" ]]; then
    N_TX=$(tail -n +2 "${OUTDIR}/abundance.tsv" | wc -l)
    echo ""
    echo "[$(date)] Quantification complete in ${ELAPSED}s"
    echo "  Transcripts quantified: ${N_TX}"
    echo "  Output: ${OUTDIR}/abundance.tsv"

    # Extract total reads from run_info.json
    if [[ -f "${OUTDIR}/run_info.json" ]]; then
        N_PROC=$(grep -o '"n_processed": [0-9]*' "${OUTDIR}/run_info.json" | grep -o '[0-9]*')
        N_PSEUDO=$(grep -o '"n_pseudoaligned": [0-9]*' "${OUTDIR}/run_info.json" | grep -o '[0-9]*')
        PCT=$(grep -o '"p_pseudoaligned": [0-9.]*' "${OUTDIR}/run_info.json" | grep -o '[0-9.]*')
        echo "  Reads processed:       ${N_PROC}"
        echo "  Reads pseudoaligned:   ${N_PSEUDO} (${PCT}%)"
    fi
else
    echo "ERROR: Kallisto output not found. Check logs."
    exit 1
fi
