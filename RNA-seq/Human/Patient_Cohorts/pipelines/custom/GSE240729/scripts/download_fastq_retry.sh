#!/bin/bash
#SBATCH --job-name=liver_retry_GSE240729
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/logs/download/retry_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/logs/download/retry_%A_%a.err
#SBATCH --partition=cpu,io
#SBATCH --time=24:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16

set -euo pipefail

MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/scripts"
FASTQ_DIR="/scratch/jameslee/MASLD_downloads/GSE240729"
FINAL_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE240729/fastq"
mkdir -p "${FASTQ_DIR}" "${FINAL_DIR}"

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
ACCESSION=$(sed -n "${TASK_ID}p" "${SCRIPT_DIR}/retry_accession_list.txt")

if [[ -z "$ACCESSION" ]]; then
    echo "No accession for task ${TASK_ID}"
    exit 0
fi

echo "=== Retry download: ${ACCESSION} (task ${TASK_ID}) ==="
echo "Start time: $(date)"

# Download
echo "Starting fasterq-dump for ${ACCESSION}..."
fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 16 --force --progress "$ACCESSION"
echo "fasterq-dump complete: $(date)"

# Compress
echo "Compressing ${ACCESSION}..."
gzip -f "${FASTQ_DIR}/${ACCESSION}"*.fastq
echo "Compression complete: $(date)"

# Verify integrity
echo "Verifying gzip integrity..."
for f in "${FASTQ_DIR}/${ACCESSION}"*.fastq.gz; do
    if ! gzip -t "$f"; then
        echo "CORRUPT: $f"
        exit 1
    fi
    echo "  OK: $(basename $f) ($(du -h $f | cut -f1))"
done

# Move to final destination
echo "Moving to final directory..."
mv "${FASTQ_DIR}/${ACCESSION}"*.fastq.gz "${FINAL_DIR}/"
echo "Done: ${ACCESSION} at $(date)"
