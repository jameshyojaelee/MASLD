#!/bin/bash
#SBATCH --job-name=validate_hra007511
#SBATCH --output=scripts/data_retrieval/logs/%x_%A.out
#SBATCH --error=scripts/data_retrieval/logs/%x_%A.err
#SBATCH --partition=cpu,io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00

OUT_DIR="data/Spatial/HRA007511/fastq"
FAILED_LOG="${OUT_DIR}/failed_fastqs.log"
> "$FAILED_LOG"

echo "Validating fastq.gz files in $OUT_DIR"

# Ensure directory exists and has files
if [ ! -d "$OUT_DIR" ] || [ -z "$(ls -A "$OUT_DIR")" ]; then
    echo "Directory $OUT_DIR is empty or missing! Marking for rescue..."
    echo "missing_directory" >> "$FAILED_LOG"
else
    # Find all gzipped fastq files (BGI formats them inside subdirectories usually)
    find "$OUT_DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | while read -r fq; do
        if ! gzip -t "$fq" 2>/dev/null; then
            echo "Corrupted file detected: $fq"
            # Delete corrupted partial download so wget -c resets cleanly
            rm -f "$fq"
            echo "$fq" >> "$FAILED_LOG"
        fi
    done
fi

if [ -s "$FAILED_LOG" ]; then
    NUM_FAILED=$(wc -l < "$FAILED_LOG")
    echo "------------------------------------------------------"
    echo "Found $NUM_FAILED corrupted files or missing root. Resubmitting wget download..."
    echo "Failed check lines: $(cat "$FAILED_LOG" | head -n 5)..."
    echo "------------------------------------------------------"
    
    # Resubmit the main download job
    RESCUE_JOB_ID=$(sbatch --parsable scripts/data_retrieval/download_hra007511.sh)
    echo "Submitted rescue download job: $RESCUE_JOB_ID"
    
    # Recursively submit this validation script to check the rescued files
    sbatch --dependency=afterany:$RESCUE_JOB_ID scripts/data_retrieval/validate_hra007511.sh
else
    echo "------------------------------------------------------"
    echo "SUCCESS: All HRA007511 spatial data files validated successfully (gzip structure intact)."
    echo "------------------------------------------------------"
    rm -f "$FAILED_LOG"
fi
