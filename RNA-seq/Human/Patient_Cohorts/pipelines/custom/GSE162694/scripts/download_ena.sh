#!/bin/bash
#SBATCH --job-name=liver_dl_ena_GSE162694
#SBATCH --output=download_%x_%A_%a.out
#SBATCH --error=download_%x_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=2

set -euo pipefail

FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE162694/fastq"
mkdir -p "${FASTQ_DIR}"

CHUNK_SIZE=5
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
START_LINE=$(( ($TASK_ID - 1) * $CHUNK_SIZE + 1 ))
END_LINE=$(( $TASK_ID * $CHUNK_SIZE ))

echo "Job Array ID: $TASK_ID"
echo "Processing lines $START_LINE to $END_LINE from ena_download_list.tsv"

sed -n "${START_LINE},${END_LINE}p" ena_download_list.tsv | while IFS=$'\t' read -r ACCESSION URL; do
    if [[ -z "$ACCESSION" || -z "$URL" ]]; then continue; fi

    OUTFILE="${FASTQ_DIR}/${ACCESSION}.fastq.gz"

    # Skip if already downloaded
    if [[ -f "$OUTFILE" ]]; then
        echo "Already exists: $ACCESSION — skipping"
        continue
    fi

    echo "Downloading $ACCESSION from ENA..."
    if ! wget -q --tries=3 --timeout=300 -O "$OUTFILE" "$URL"; then
        echo "Error: wget failed for $ACCESSION, removing partial file"
        rm -f "$OUTFILE"
        exit 1
    fi

    # Verify gzip integrity
    if ! gzip -t "$OUTFILE"; then
        echo "Error: corrupt gzip for $ACCESSION, removing"
        rm -f "$OUTFILE"
        exit 1
    fi

    echo "OK: $ACCESSION ($(stat --printf='%s' "$OUTFILE") bytes)"
done

echo "Chunk complete."
