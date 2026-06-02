#!/bin/bash
#SBATCH --job-name=liver_dl2_GSE174478
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

FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE174478/fastq"
mkdir -p "${FASTQ_DIR}"

CHUNK_SIZE=5
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
START_LINE=$(( ($TASK_ID - 1) * $CHUNK_SIZE + 1 ))
END_LINE=$(( $TASK_ID * $CHUNK_SIZE ))

echo "Job Array ID: $TASK_ID"
echo "Processing lines $START_LINE to $END_LINE from accession_list_retry2.txt"

SCRATCH_TEMP="/scratch/fasterq_${SLURM_JOB_ID}_${TASK_ID}"
mkdir -p "${SCRATCH_TEMP}"
trap "rm -rf ${SCRATCH_TEMP}" EXIT

sed -n "${START_LINE},${END_LINE}p" accession_list_retry2.txt | while read -r ACCESSION; do
    if [[ -z "$ACCESSION" ]]; then continue; fi

    # Skip if already downloaded
    if ls "${FASTQ_DIR}/${ACCESSION}"*.fastq.gz 1>/dev/null 2>&1; then
        echo "Already exists: $ACCESSION — skipping"
        continue
    fi

    echo "Prefetching $ACCESSION..."
    prefetch --max-size 50G --output-directory "${SCRATCH_TEMP}" "$ACCESSION"

    echo "Running fasterq-dump on local SRA for $ACCESSION..."
    if ! fasterq-dump --outdir "${FASTQ_DIR}" --temp "${SCRATCH_TEMP}" --split-files --threads 8 --force --progress "${SCRATCH_TEMP}/${ACCESSION}/${ACCESSION}.sra"; then
        echo "Error: fasterq-dump failed for $ACCESSION"
        exit 1
    fi

    echo "Compressing ${ACCESSION}..."
    gzip -f "${FASTQ_DIR}/${ACCESSION}"*.fastq

    rm -rf "${SCRATCH_TEMP}/${ACCESSION}"
done

wait
echo "Chunk download and compression complete."
