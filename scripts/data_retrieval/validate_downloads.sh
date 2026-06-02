#!/bin/bash
#SBATCH --job-name=validate_fastq
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4:00:00

# Usage: sbatch --array=1-N validate_downloads.sh <SRR_LIST_FILE> <OUTPUT_DIR>

SRR_LIST=$1
OUT_DIR=$2

if [ -z "$SRR_LIST" ] || [ -z "$OUT_DIR" ]; then
    echo "Error: Must provide SRR_LIST and OUT_DIR"
    exit 1
fi

SRR_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")

if [ -z "$SRR_ID" ]; then
    echo "Error: Could not determine SRR ID for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Validating $SRR_ID..."
TARGET_DIR="$OUT_DIR/$SRR_ID"

if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Directory $TARGET_DIR does not exist."
    exit 1
fi

# Check gzip integrity
count=0
for fq in "$TARGET_DIR"/*.fastq.gz; do
    if [ -f "$fq" ]; then
        echo "Testing integrity of $fq..."
        count=$((count+1))
        gzip -t "$fq"
        if [ $? -ne 0 ]; then
            echo "ERROR: Corrupted gzip file detected: $fq"
            exit 1
        fi
    fi
done

if [ "$count" -eq 0 ]; then
    echo "ERROR: No fastq.gz files found in $TARGET_DIR"
    exit 1
fi

echo "SUCCESS: All files in $TARGET_DIR validated successfully."
