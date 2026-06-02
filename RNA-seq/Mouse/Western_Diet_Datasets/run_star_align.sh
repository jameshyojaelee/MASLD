#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=logs/star_%A_%a.out
#SBATCH --error=logs/star_%A_%a.err
# Submit per dataset: sbatch --job-name=STAR_GSE220575 --array=1-17%25 run_star_align.sh GSE220575

set -eo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

DATASET_ID="${1:?Usage: sbatch --array=1-N%25 run_star_align.sh DATASET_ID}"
PROJECT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WD_DIR="${PROJECT}/RNA-seq/Mouse/Western_Diet_Datasets"
STAR_INDEX="${PROJECT}/RNA-seq/Mouse/Public_Diet_Models/reference/star_vm38"
SAMPLE_LIST="${WD_DIR}/${DATASET_ID}/metadata/sample_list.txt"

TASK_ID=${SLURM_ARRAY_TASK_ID}
LINE=$(sed -n "${TASK_ID}p" "${SAMPLE_LIST}")
SRR=$(echo "${LINE}" | cut -f1)
NCOLS=$(echo "${LINE}" | awk -F'\t' '{print NF}')

resolve_path() {
    local p="$1"
    if [[ "$p" = /* ]]; then echo "$p"; else echo "${PROJECT}/${p}"; fi
}

R1=$(resolve_path "$(echo "${LINE}" | cut -f2)")

OUTDIR="${WD_DIR}/${DATASET_ID}/alignments/star/${SRR}"
mkdir -p "${OUTDIR}"

echo "=== STAR alignment: ${SRR} (dataset: ${DATASET_ID}) ==="
echo "Task ${TASK_ID}, PE columns: ${NCOLS}"

if [[ "${NCOLS}" -ge 3 ]]; then
    R2=$(resolve_path "$(echo "${LINE}" | cut -f3)")
    echo "PE: R1=${R1} R2=${R2}"
    STAR --runThreadN ${SLURM_CPUS_PER_TASK} \
         --genomeDir "${STAR_INDEX}" \
         --readFilesIn "${R1}" "${R2}" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${OUTDIR}/${SRR}." \
         --outSAMtype BAM SortedByCoordinate \
         --twopassMode Basic \
         --quantMode GeneCounts \
         --limitBAMsortRAM 15000000000
else
    echo "SE: R1=${R1}"
    STAR --runThreadN ${SLURM_CPUS_PER_TASK} \
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
    echo "SUCCESS: BAM created"
    UNIQ=$(grep "Uniquely mapped reads %" "${OUTDIR}/${SRR}.Log.final.out" | awk '{print $NF}')
    echo "Uniquely mapped: ${UNIQ}"
else
    echo "FAILED: no BAM output" >&2
    exit 1
fi
