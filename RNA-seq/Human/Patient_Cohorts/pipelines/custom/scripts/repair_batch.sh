#!/bin/bash
#SBATCH --job-name=liver_repair_batch
#SBATCH --output=logs/repair_batch_%j.out
#SBATCH --error=logs/repair_batch_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-20

set -euo pipefail

# USAGE: sbatch repair_batch.sh <BAD_FILES_LIST>

BAD_FILES_LIST=$1

if [[ ! -f "$BAD_FILES_LIST" ]]; then
    echo "Error: List $BAD_FILES_LIST not found."
    exit 1
fi

# Setup Environment
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

# Process Chunk
CHUNK_SIZE=5
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
START_LINE=$(( ($TASK_ID - 1) * $CHUNK_SIZE + 1 ))
END_LINE=$(( $TASK_ID * $CHUNK_SIZE ))

echo "Processing lines $START_LINE to $END_LINE from $BAD_FILES_LIST"

sed -n "${START_LINE},${END_LINE}p" "$BAD_FILES_LIST" | while read -r FILEPATH; do
    if [[ -z "$FILEPATH" ]]; then continue; fi
    
    # Extract Accession (SRRXXXXXX) from filename
    FILENAME=$(basename "$FILEPATH")
    ACCESSION=$(echo "$FILENAME" | grep -o "SRR[0-9]*")
    DIRNAME=$(dirname "$FILEPATH")
    
    echo "Repairing $ACCESSION in $DIRNAME..."
    
    # Remove corrupt file
    rm -f "$FILEPATH"
    
    # Re-download
    fasterq-dump --outdir "$DIRNAME" --split-files --threads 4 --force --progress "$ACCESSION"
    
    # Compress
    gzip -f "$DIRNAME"/${ACCESSION}*.fastq
    
    echo "Verified: $FILEPATH repaired."
done
