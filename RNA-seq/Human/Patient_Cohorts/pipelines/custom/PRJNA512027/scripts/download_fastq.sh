#!/bin/bash
#SBATCH --job-name=liver_dl_PRJNA512027
#SBATCH --output=download_%x_%A_%a.out
#SBATCH --error=download_%x_%A_%a.err
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8

#SBATCH --array=1-39

set -euo pipefail

# Activate environment with sra-tools
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

PROJECT_ROOT=$(pwd)
METADATA_DIR="../metadata"
RESULTS_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/PRJNA512027"

FASTQ_DIR="${RESULTS_DIR}/fastq"
mkdir -p "${FASTQ_DIR}"

INPUT_FILE=${1:-accession_list.txt}

if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "Error: ${INPUT_FILE} not found in $(pwd)"
    exit 1
fi

# Array Job Logic
CHUNK_SIZE=5
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}

START_LINE=$(( ($TASK_ID - 1) * $CHUNK_SIZE + 1 ))
END_LINE=$(( $TASK_ID * $CHUNK_SIZE ))

echo "Job Array ID: $TASK_ID"
echo "Processing lines $START_LINE to $END_LINE from ${INPUT_FILE}"

sed -n "${START_LINE},${END_LINE}p" "${INPUT_FILE}" | while read -r ACCESSION; do
    if [[ -z "$ACCESSION" ]]; then continue; fi
    
    # Check if files already exist (simple check for _1.fastq.gz and _2.fastq.gz)
    if [[ -f "${FASTQ_DIR}/${ACCESSION}_1.fastq.gz" && -f "${FASTQ_DIR}/${ACCESSION}_2.fastq.gz" ]]; then
        echo "Skipping $ACCESSION (Files already exist)"
        continue
    fi

    echo "Starting download for $ACCESSION..."
    # Removed --force to be safer, though logic above handles skipping
    if ! fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 8 --progress "$ACCESSION"; then
        echo "Error: fasterq-dump failed for $ACCESSION"
        # Don't exit immediately, try next sample in chunk
        continue
    fi
    
    echo "Compressing ${ACCESSION}..."
    gzip -f "${FASTQ_DIR}/${ACCESSION}"*.fastq

done

echo "Chunk download and compression complete."
