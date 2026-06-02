#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=logs/star_bigmem_%j.out
#SBATCH --error=logs/star_bigmem_%j.err
# Runs multiple STAR alignments in parallel on one bigmem node.
# 64 CPUs → 4 concurrent STAR jobs × 16 threads each.
# Usage: sbatch --job-name=STAR_BM_GSE246328 run_star_bigmem.sh GSE246328

set -eo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

DATASET_ID="${1:?Usage: sbatch run_star_bigmem.sh DATASET_ID}"
PROJECT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WD_DIR="${PROJECT}/RNA-seq/Mouse/Western_Diet_Datasets"
STAR_INDEX="${PROJECT}/RNA-seq/Mouse/Public_Diet_Models/reference/star_vm38"
SAMPLE_LIST="${WD_DIR}/${DATASET_ID}/metadata/sample_list.txt"
NSAMPLES=$(wc -l < "${SAMPLE_LIST}")
TOTAL_CPUS=${SLURM_CPUS_PER_TASK}
THREADS_PER=$(( TOTAL_CPUS < 16 ? TOTAL_CPUS : 16 ))
PARALLEL_JOBS=$(( TOTAL_CPUS / THREADS_PER ))
if [[ ${PARALLEL_JOBS} -lt 1 ]]; then PARALLEL_JOBS=1; fi

resolve_path() {
    local p="$1"
    if [[ "$p" = /* ]]; then echo "$p"; else echo "${PROJECT}/${p}"; fi
}

align_sample() {
    local TASK_ID=$1
    local LINE=$(sed -n "${TASK_ID}p" "${SAMPLE_LIST}")
    local SRR=$(echo "${LINE}" | cut -f1)
    local NCOLS=$(echo "${LINE}" | awk -F'\t' '{print NF}')
    local R1=$(resolve_path "$(echo "${LINE}" | cut -f2)")
    local OUTDIR="${WD_DIR}/${DATASET_ID}/alignments/star/${SRR}"

    if [[ -f "${OUTDIR}/${SRR}.Aligned.sortedByCoord.out.bam" ]]; then
        echo "[${TASK_ID}/${NSAMPLES}] ${SRR}: BAM exists, skipping"
        return 0
    fi

    mkdir -p "${OUTDIR}"

    if [[ "${NCOLS}" -ge 3 ]]; then
        local R2=$(resolve_path "$(echo "${LINE}" | cut -f3)")
        STAR --runThreadN ${THREADS_PER} \
             --genomeDir "${STAR_INDEX}" \
             --readFilesIn "${R1}" "${R2}" \
             --readFilesCommand zcat \
             --outFileNamePrefix "${OUTDIR}/${SRR}." \
             --outSAMtype BAM SortedByCoordinate \
             --twopassMode Basic \
             --quantMode GeneCounts \
             --limitBAMsortRAM 15000000000
    else
        STAR --runThreadN ${THREADS_PER} \
             --genomeDir "${STAR_INDEX}" \
             --readFilesIn "${R1}" \
             --readFilesCommand zcat \
             --outFileNamePrefix "${OUTDIR}/${SRR}." \
             --outSAMtype BAM SortedByCoordinate \
             --twopassMode Basic \
             --quantMode GeneCounts \
             --limitBAMsortRAM 15000000000
    fi

    if [[ -f "${OUTDIR}/${SRR}.Aligned.sortedByCoord.out.bam" ]]; then
        local UNIQ=$(grep "Uniquely mapped reads %" "${OUTDIR}/${SRR}.Log.final.out" 2>/dev/null | awk '{print $NF}')
        echo "[${TASK_ID}/${NSAMPLES}] ${SRR}: OK (${UNIQ})"
    else
        echo "[${TASK_ID}/${NSAMPLES}] ${SRR}: FAILED" >&2
        return 1
    fi
}
export -f align_sample resolve_path
export SAMPLE_LIST NSAMPLES WD_DIR DATASET_ID STAR_INDEX THREADS_PER PROJECT

echo "=== STAR bigmem: ${DATASET_ID} (${NSAMPLES} samples, ${PARALLEL_JOBS} parallel × ${THREADS_PER} threads) ==="

# Run in batches of PARALLEL_JOBS using bash background jobs
for BATCH_START in $(seq 1 ${PARALLEL_JOBS} ${NSAMPLES}); do
    BATCH_END=$(( BATCH_START + PARALLEL_JOBS - 1 ))
    if [[ ${BATCH_END} -gt ${NSAMPLES} ]]; then BATCH_END=${NSAMPLES}; fi

    PIDS=()
    for TASK_ID in $(seq ${BATCH_START} ${BATCH_END}); do
        align_sample ${TASK_ID} &
        PIDS+=($!)
    done

    for PID in "${PIDS[@]}"; do
        wait ${PID} || true
    done
    echo "Batch ${BATCH_START}-${BATCH_END} complete"
done

DONE=$(find "${WD_DIR}/${DATASET_ID}/alignments/star" -name "*.Aligned.sortedByCoord.out.bam" 2>/dev/null | wc -l)
echo "=== ${DATASET_ID} DONE: ${DONE}/${NSAMPLES} BAMs ==="
