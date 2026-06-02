#!/bin/bash
#SBATCH --job-name=liver_rescue_dl_array
#SBATCH --output=logs/rescue_array_%A_%a.out
#SBATCH --error=logs/rescue_array_%A_%a.err
#SBATCH --partition=io
#SBATCH --array=0-5
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
FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE213621/fastq"
mkdir -p "$FASTQ_DIR"

# List of corrupted accessions
ACCESSIONS=(
"SRR21622948"
"SRR21622952"
"SRR21623064"
"SRR21623089"
"SRR21623146"
"SRR21623152"
)

# Select accession based on Task ID
ACCESSION=${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}

echo "=================================================="
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Rescuing Accession: $ACCESSION"

# 1. Remove potential corrupted/partial files
echo "Cleaning up old files for $ACCESSION..."
rm -f "${FASTQ_DIR}/${ACCESSION}"*.fastq.gz
rm -f "${FASTQ_DIR}/${ACCESSION}"*.fastq

# 2. Re-download (Multithreaded)
echo "Downloading $ACCESSION with fasterq-dump..."
# fasterq-dump --threads 8 uses ~8 cores
fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 8 --force --progress "$ACCESSION"

# 3. Compress (Parallel gzip if available, otherwise gzip)
echo "Compressing $ACCESSION..."
if command -v pigz &> /dev/null; then
    echo "Using pigz for parallel compression..."
    # pigz uses all available cores by default, limiting to 8 to match allocation
    pigz -p 8 "${FASTQ_DIR}/${ACCESSION}"*.fastq
else
    echo "pigz not found, using standard gzip..."
    gzip "${FASTQ_DIR}/${ACCESSION}"*.fastq
fi

echo "Finished $ACCESSION"
