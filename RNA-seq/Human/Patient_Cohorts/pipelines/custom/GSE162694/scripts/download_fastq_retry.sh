#!/bin/bash
#SBATCH --job-name=liver_retry_GSE162694
#SBATCH --output=download_%x_%A_%a.out
#SBATCH --error=download_%x_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8

set -euo pipefail

MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE162694/fastq"
mkdir -p "${FASTQ_DIR}"

CHUNK_SIZE=5
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
START_LINE=$(( ($TASK_ID - 1) * $CHUNK_SIZE + 1 ))
END_LINE=$(( $TASK_ID * $CHUNK_SIZE ))

echo "Job Array ID: $TASK_ID"
echo "Processing lines $START_LINE to $END_LINE from accession_list_retry.txt"

SCRATCH_TEMP="/scratch/fasterq_${SLURM_JOB_ID}_${TASK_ID}"
mkdir -p "${SCRATCH_TEMP}"
trap "rm -rf ${SCRATCH_TEMP}" EXIT

sed -n "${START_LINE},${END_LINE}p" accession_list_retry.txt | while read -r ACCESSION; do
    if [[ -z "$ACCESSION" ]]; then continue; fi
    echo "Starting download for $ACCESSION..."
    if ! fasterq-dump --outdir "${FASTQ_DIR}" --temp "${SCRATCH_TEMP}" --split-files --threads 8 --force --progress "$ACCESSION"; then
        echo "Error: fasterq-dump failed for $ACCESSION"
        exit 1
    fi
    echo "Compressing ${ACCESSION}..."
    gzip -f "${FASTQ_DIR}/${ACCESSION}"*.fastq
done

wait
echo "Chunk download and compression complete."
