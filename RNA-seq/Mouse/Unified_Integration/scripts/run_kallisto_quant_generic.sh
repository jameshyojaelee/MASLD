#!/bin/bash
#SBATCH --job-name=kallisto_quant
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=logs/kallisto_quant_%A_%a.out
#SBATCH --error=logs/kallisto_quant_%A_%a.err

# =============================================================================
# Generic Kallisto quantification — works for any mouse dataset
# =============================================================================
# Usage:
#   # Paired-end:
#   sbatch --array=1-N run_kallisto_quant_generic.sh /path/to/sample_list.txt /path/to/outdir
#
#   # Single-end (fragment length + SD auto-added if missing):
#   sbatch --array=1-N run_kallisto_quant_generic.sh /path/to/sample_list.txt /path/to/outdir --single -l 200 -s 30
#
# Input sample_list.txt — tab-separated, no header:
#   PE: sample_id \t R1_path \t R2_path
#   SE: sample_id \t fastq_path
#
# Output per sample: <outdir>/<sample_id>/{abundance.h5, abundance.tsv, run_info.json}
# =============================================================================

set -euo pipefail

module load kallisto/0.51.1

# ---- Configuration ----
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INDEX="${BASE}/data/reference/kallisto/gencode.vM38.kallisto.idx"

# ---- Parse arguments ----
SAMPLE_LIST="${1:?ERROR: arg1 = sample_list.txt path}"
OUTBASE="${2:?ERROR: arg2 = output base directory}"
shift 2
EXTRA_FLAGS=("$@")

# ---- Validate ----
if [[ ! -f "${INDEX}" ]]; then
    echo "ERROR: Kallisto index not found: ${INDEX}" ; exit 1
fi
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: Sample list not found: ${SAMPLE_LIST}" ; exit 1
fi

# ---- Extract sample info for this array task ----
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")
if [[ -z "${LINE}" ]]; then
    echo "ERROR: No entry for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" ; exit 1
fi

NCOLS=$(echo "${LINE}" | awk -F'\t' '{print NF}')
SRR=$(echo "${LINE}" | cut -f1)
OUTDIR="${OUTBASE}/${SRR}"
mkdir -p "${OUTDIR}"

echo "============================================="
echo "  Kallisto Quantification"
echo "============================================="
echo "  Sample:   ${SRR} (task ${SLURM_ARRAY_TASK_ID})"
echo "  Index:    ${INDEX}"
echo "  Threads:  ${SLURM_CPUS_PER_TASK}"
echo "  Output:   ${OUTDIR}"

START_TIME=$SECONDS

if [[ ${NCOLS} -ge 3 ]]; then
    R1=$(echo "${LINE}" | cut -f2)
    R2=$(echo "${LINE}" | cut -f3)
    echo "  Mode:     paired-end"
    echo "  R1:       ${R1}"
    echo "  R2:       ${R2}"
    echo "============================================="
    for f in "${R1}" "${R2}"; do
        [[ ! -f "${f}" ]] && echo "ERROR: FASTQ not found: ${f}" && exit 1
    done
    kallisto quant -i "${INDEX}" -o "${OUTDIR}" -t "${SLURM_CPUS_PER_TASK}" \
        "${EXTRA_FLAGS[@]+"${EXTRA_FLAGS[@]}"}" "${R1}" "${R2}"

elif [[ ${NCOLS} -eq 2 ]]; then
    FQ=$(echo "${LINE}" | cut -f2)
    echo "  Mode:     single-end"
    echo "  FASTQ:    ${FQ}"
    echo "============================================="
    [[ ! -f "${FQ}" ]] && echo "ERROR: FASTQ not found: ${FQ}" && exit 1
    HAS_SINGLE=0
    for flag in "${EXTRA_FLAGS[@]+"${EXTRA_FLAGS[@]}"}"; do
        [[ "${flag}" == "--single" ]] && HAS_SINGLE=1
    done
    if [[ ${HAS_SINGLE} -eq 0 ]]; then
        echo "  Auto-adding --single -l 200 -s 30"
        EXTRA_FLAGS+=(--single -l 200 -s 30)
    fi
    kallisto quant -i "${INDEX}" -o "${OUTDIR}" -t "${SLURM_CPUS_PER_TASK}" \
        "${EXTRA_FLAGS[@]}" "${FQ}"
else
    echo "ERROR: Unexpected column count (${NCOLS}) in line ${SLURM_ARRAY_TASK_ID}" ; exit 1
fi

ELAPSED=$(( SECONDS - START_TIME ))

if [[ -f "${OUTDIR}/abundance.tsv" ]]; then
    N_TX=$(tail -n +2 "${OUTDIR}/abundance.tsv" | wc -l)
    echo ""
    echo "[$(date)] Done in ${ELAPSED}s — ${N_TX} transcripts quantified"
    if [[ -f "${OUTDIR}/run_info.json" ]]; then
        N_PROC=$(grep -o '"n_processed": [0-9]*' "${OUTDIR}/run_info.json" | grep -o '[0-9]*')
        N_PSEUDO=$(grep -o '"n_pseudoaligned": [0-9]*' "${OUTDIR}/run_info.json" | grep -o '[0-9]*')
        PCT=$(grep -o '"p_pseudoaligned": [0-9.]*' "${OUTDIR}/run_info.json" | grep -o '[0-9.]*')
        echo "  Reads processed:     ${N_PROC}"
        echo "  Pseudoaligned:       ${N_PSEUDO} (${PCT}%)"
    fi
else
    echo "ERROR: Kallisto output not found." ; exit 1
fi
