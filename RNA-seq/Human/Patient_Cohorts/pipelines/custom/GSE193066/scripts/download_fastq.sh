#!/bin/bash
#SBATCH --job-name=liver_dl_GSE193066
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066/logs/download/download_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066/logs/download/download_%A_%a.err
#SBATCH --partition=cpu,io
#SBATCH --time=12:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8

set -euo pipefail

MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

PROJECT_ROOT=$(pwd)
FASTQ_DIR="/scratch/jameslee/MASLD_downloads/GSE193066"
mkdir -p "${FASTQ_DIR}"
cd "${FASTQ_DIR}"

if [[ ! -f "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066/scripts/accession_list.txt" ]]; then
    echo "Error: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066/scripts/accession_list.txt not found"
    exit 1
fi

CHUNK_SIZE=1
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}

START_LINE=$(( ($TASK_ID - 1) * $CHUNK_SIZE + 1 ))
END_LINE=$(( $TASK_ID * $CHUNK_SIZE ))

sed -n "${START_LINE},${END_LINE}p" /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066/scripts/accession_list.txt | while read -r ACCESSION; do
    if [[ -z "$ACCESSION" ]]; then continue; fi
    echo "Starting download for $ACCESSION..."
    if ! fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 8 --force --progress "$ACCESSION"; then
        echo "Error: fasterq-dump failed for $ACCESSION"
        exit 1
    fi
    echo "Compressing ${ACCESSION}..."
    gzip -f "${FASTQ_DIR}/${ACCESSION}"*.fastq
    echo "Moving ${ACCESSION} to project directory..."
    mv "${FASTQ_DIR}/${ACCESSION}"*.fastq.gz /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE193066/fastq/
done
wait
echo "Chunk download and compression complete."
