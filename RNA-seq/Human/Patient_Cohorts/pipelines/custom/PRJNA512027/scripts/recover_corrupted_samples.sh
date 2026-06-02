#!/bin/bash
#SBATCH --job-name=liver_recover_prjna
#SBATCH --output=logs/recover_prjna_%j.out
#SBATCH --error=logs/recover_prjna_%j.err
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8

set -euo pipefail

# Setup
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

PROJECT_ROOT=$(pwd)
METADATA_DIR="pipelines/custom/PRJNA512027/metadata"
RESULTS_DIR="data/raw/PRJNA512027"
FASTQ_DIR="${RESULTS_DIR}/fastq"
INPUT_FILE="${METADATA_DIR}/recovery_accessions.txt"

if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "Error: ${INPUT_FILE} not found."
    exit 1
fi

echo "Starting recovery for samples in ${INPUT_FILE}..."

while read -r ACCESSION; do
    if [[ -z "$ACCESSION" ]]; then continue; fi

    echo "------------------------------------------------"
    echo "Recovering $ACCESSION"
    
    # 1. Delete Corrupted Files
    echo "Deleting corrupted files for $ACCESSION..."
    rm -f "${FASTQ_DIR}/${ACCESSION}"*.fastq.gz "${FASTQ_DIR}/${ACCESSION}"*.fastq
    
    # 2. Re-download
    echo "Downloading $ACCESSION..."
    if ! fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 8 --progress "$ACCESSION"; then
        echo "Error: fasterq-dump failed for $ACCESSION"
        continue
    fi
    
    # 3. Compress
    echo "Compressing $ACCESSION..."
    gzip -f "${FASTQ_DIR}/${ACCESSION}"*.fastq
    
    echo "Success: $ACCESSION recovered."

done < "${INPUT_FILE}"

echo "Recovery job complete."
