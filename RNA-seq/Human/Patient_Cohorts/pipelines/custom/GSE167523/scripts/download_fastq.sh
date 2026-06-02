#!/bin/bash
#SBATCH --job-name=liver_dl_GSE167523
#SBATCH --output=download_%x_%A_%a.out
#SBATCH --error=download_%x_%A_%a.err
#SBATCH --partition=io
#SBATCH --time=12:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8


set -euo pipefail

# Activate environment with sra-tools
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

PROJECT_ROOT=$(pwd)
METADATA_DIR="../metadata"
RESULTS_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE167523"

FASTQ_DIR="${RESULTS_DIR}/fastq"
mkdir -p "${FASTQ_DIR}"

if [[ ! -f "accession_list.txt" ]]; then
    echo "Error: accession_list.txt not found in $(pwd)"
    exit 1
fi

# Array Job Logic
CHUNK_SIZE=5
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}

START_LINE=$(( ($TASK_ID - 1) * $CHUNK_SIZE + 1 ))
END_LINE=$(( $TASK_ID * $CHUNK_SIZE ))

echo "Job Array ID: $TASK_ID"
echo "Processing lines $START_LINE to $END_LINE from accession_list.txt"

sed -n "${START_LINE},${END_LINE}p" accession_list.txt | while read -r ACCESSION; do
    if [[ -z "$ACCESSION" ]]; then continue; fi
    
    echo "Starting download for $ACCESSION..."
    if ! fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 8 --force --progress "$ACCESSION"; then
        echo "Error: fasterq-dump failed for $ACCESSION"
        exit 1
    fi
    
    echo "Compressing ${ACCESSION}..."
    gzip -f "${FASTQ_DIR}/${ACCESSION}"*.fastq

    
done

echo "Chunk download and compression complete."
