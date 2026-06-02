#!/bin/bash
#SBATCH --job-name=liver_dl_missing4_GSE126848
#SBATCH --output=download_missing4_%j.out
#SBATCH --error=download_missing4_%j.err
#SBATCH --partition=io
#SBATCH --time=06:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8

set -euo pipefail

MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE126848/fastq"
mkdir -p "${FASTQ_DIR}"

ACCESSION_FILE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE126848/scripts/missing_4.txt"

if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "Error: missing_4.txt not found at $ACCESSION_FILE"
    exit 1
fi

while read -r ACCESSION; do
    if [[ -z "$ACCESSION" ]]; then continue; fi

    echo "Downloading $ACCESSION..."
    if ! fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 8 --force --progress "$ACCESSION"; then
        echo "Warning: fasterq-dump failed for $ACCESSION — skipping (possible SRA corruption)"
        continue
    fi

    echo "Compressing ${ACCESSION}..."
    gzip -f "${FASTQ_DIR}/${ACCESSION}"*.fastq

    echo "Done: $ACCESSION"
done < "$ACCESSION_FILE"

wait
echo "All 4 missing samples downloaded and compressed."
