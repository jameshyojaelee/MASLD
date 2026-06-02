#!/bin/bash
#SBATCH --job-name=liver_rescue_dl_GSE213621
#SBATCH --output=logs/rescue_download_%j.out
#SBATCH --error=logs/rescue_download_%j.err
#SBATCH --partition=io
#SBATCH --time=12:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8

set -euo pipefail

# Setup Environment
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

# Define paths
PROJECT_ROOT=$(dirname $(pwd)) # Assumes script is in scripts/
FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE213621/fastq"

# List of corrupted accessions
ACCESSIONS=(
"SRR21622948"
"SRR21622952"
"SRR21623064"
"SRR21623089"
"SRR21623146"
"SRR21623152"
)

mkdir -p "$FASTQ_DIR"

for ACCESSION in "${ACCESSIONS[@]}"; do
    echo "=================================================="
    echo "Rescuing $ACCESSION"
    
    # 1. Remove corrupted files
    echo "Removing potential corrupted files for $ACCESSION..."
    rm -f "${FASTQ_DIR}/${ACCESSION}"*.fastq.gz
    rm -f "${FASTQ_DIR}/${ACCESSION}"*.fastq

    # 2. Re-download
    echo "Downloading $ACCESSION..."
    fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 8 --force --progress "$ACCESSION"

    # 3. Compress
    echo "Compressing $ACCESSION..."
    gzip "${FASTQ_DIR}/${ACCESSION}"*.fastq
    
    echo "Finished $ACCESSION"
done

echo "All rescue downloads complete."
