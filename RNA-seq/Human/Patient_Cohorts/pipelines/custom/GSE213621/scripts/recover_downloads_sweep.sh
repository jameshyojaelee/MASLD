#!/bin/bash
#SBATCH --job-name=liver_recover_sweep_GSE213621
#SBATCH --output=logs/recovery_%A_%a.out
#SBATCH --error=logs/recovery_%A_%a.err
#SBATCH --partition=io
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

# USAGE: sbatch --array=1-29 scripts/recover_downloads.sh

set -euo pipefail

export TMPDIR=/scratch
export TEMP=/scratch
export TMP=/scratch

# Activate environment with sra-tools
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

PROJECT_ROOT=$(pwd)
FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE213621/fastq"
ACCESSION_FILE="scripts/recovery_accessions_sweep.txt"

if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "Error: $ACCESSION_FILE not found"
    exit 1
fi

# Array Job Logic: One accession per task
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
ACCESSION=$(sed -n "${TASK_ID}p" "$ACCESSION_FILE")

if [[ -z "$ACCESSION" ]]; then
    echo "Error: No accession found at line $TASK_ID"
    exit 1
fi

echo "------------------------------------------------"
echo "Job Array ID: $TASK_ID"
echo "Recovering $ACCESSION..."

# Use --force to overwrite corrupted files
if ! fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 8 --force --progress "$ACCESSION"; then
    echo "Error: fasterq-dump failed for $ACCESSION"
    exit 1
fi

echo "Compressing ${ACCESSION}..."
gzip -f "${FASTQ_DIR}/${ACCESSION}"*.fastq

echo "Recovery for $ACCESSION complete."
