#!/bin/bash
#SBATCH --job-name=liver_retry_dl_GSE213621
#SBATCH --output=download_%x_%A_%a.out
#SBATCH --error=download_%x_%A_%a.err
#SBATCH --partition=cpu
#SBATCH --time=24:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8

# Retry download script — uses retry_accession_list.txt
# - Longer wall time (24h vs 12h)
# - Downloads one sample per task (no chunking) so a single failure doesn't kill others

set -euo pipefail

# Use /scratch for temp files (fasterq-dump writes large intermediates)
export TMPDIR=/scratch
export TEMP=/scratch
export TMP=/scratch

MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE213621/fastq"
mkdir -p "${FASTQ_DIR}"

ACCESSION_FILE="retry_accession_list.txt"
if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "Error: $ACCESSION_FILE not found in $(pwd)"
    exit 1
fi

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
ACCESSION=$(sed -n "${TASK_ID}p" "$ACCESSION_FILE")

if [[ -z "$ACCESSION" ]]; then
    echo "No accession at line $TASK_ID"
    exit 0
fi

# Skip if already downloaded
if [[ -f "${FASTQ_DIR}/${ACCESSION}_1.fastq.gz" ]]; then
    echo "$ACCESSION already downloaded, skipping."
    exit 0
fi

echo "Downloading $ACCESSION (task $TASK_ID)..."

if ! fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 8 --force --progress "$ACCESSION"; then
    echo "Error: fasterq-dump failed for $ACCESSION"
    exit 1
fi

echo "Compressing ${ACCESSION}..."
gzip -f "${FASTQ_DIR}/${ACCESSION}"*.fastq

echo "Done: $ACCESSION"
