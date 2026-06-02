#!/bin/bash
#SBATCH --job-name=validate_and_rescue_sra
#SBATCH --output=scripts/data_retrieval/logs/%x_%A_%a.out
#SBATCH --error=scripts/data_retrieval/logs/%x_%A_%a.err
#SBATCH --partition=cpu,io
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=12:00:00

# Usage: sbatch validate_and_rescue.sh <SRR_LIST_FILE> <OUTPUT_DIR> <DATASET_NAME> [ATTEMPT]
# The ATTEMPT counter prevents infinite recursion (max 3 retries).
#
# T3.7 safety (2026-04-22): the destructive `rm -rf "$TARGET_DIR"` on corrupted
# downloads now requires DRY_RUN=false. Default is dry-run so operators review
# what would be deleted before committing. Set DRY_RUN=false (or pass
# --execute as first argument) to actually delete corrupted directories.

# Allow --execute (or --dry-run) as the first arg to toggle behaviour.
DRY_RUN="${DRY_RUN:-true}"
if [ "${1:-}" = "--execute" ]; then
    DRY_RUN=false
    shift
elif [ "${1:-}" = "--dry-run" ]; then
    DRY_RUN=true
    shift
fi

SRR_LIST=$1
OUT_DIR=$2
DATASET_NAME=$3
ATTEMPT=${4:-1}
MAX_RETRIES=3

if [ -z "$SRR_LIST" ] || [ -z "$OUT_DIR" ] || [ -z "$DATASET_NAME" ]; then
    echo "Error: Must provide SRR_LIST, OUT_DIR, and DATASET_NAME"
    exit 1
fi

echo "=========================================="
echo "Validation attempt $ATTEMPT / $MAX_RETRIES for $DATASET_NAME"
echo "SRR list: $SRR_LIST"
echo "Output directory: $OUT_DIR"
echo "DRY_RUN mode: $DRY_RUN (set DRY_RUN=false or pass --execute to delete)"
echo "=========================================="

# Use a stable failed-list name (overwritten each attempt, not appended)
FAILED_LIST="${SRR_LIST%.txt}_rescue.txt"
> "$FAILED_LIST"

while read -r SRR_ID; do
    if [ -z "$SRR_ID" ]; then continue; fi
    TARGET_DIR="$OUT_DIR/$SRR_ID"

    FAILED=0

    if [ ! -d "$TARGET_DIR" ]; then
        echo "FAIL: Missing directory for $SRR_ID"
        FAILED=1
    else
        # Count fastq.gz files
        count=$(ls -1 "$TARGET_DIR"/*.fastq.gz 2>/dev/null | wc -l)
        if [ "$count" -eq 0 ]; then
            echo "FAIL: No fastq.gz files found for $SRR_ID"
            FAILED=1
        else
            for fq in "$TARGET_DIR"/*.fastq.gz; do
                gzip -t "$fq" 2>/dev/null
                if [ $? -ne 0 ]; then
                    echo "FAIL: Corrupted file detected: $fq"
                    FAILED=1
                fi
            done
            if [ "$FAILED" -eq 0 ]; then
                echo "OK: $SRR_ID ($count fastq.gz files)"
            fi
        fi
    fi

    if [ "$FAILED" -eq 1 ]; then
        echo "$SRR_ID" >> "$FAILED_LIST"
        # Clean up the corrupted directory so prefetch/fasterq-dump starts fresh.
        # T3.7 (2026-04-22): gated behind DRY_RUN to prevent accidental data
        # deletion. Default mode is dry-run: pass --execute (or export
        # DRY_RUN=false) to actually remove corrupted directories.
        if [ "$DRY_RUN" = "false" ]; then
            echo "  [EXECUTE] rm -rf $TARGET_DIR"
            rm -rf "$TARGET_DIR"
        else
            echo "  [DRY_RUN] would rm -rf $TARGET_DIR"
        fi
    fi
done < "$SRR_LIST"

NUM_FAILED=$(wc -l < "$FAILED_LIST")
TOTAL=$(wc -l < "$SRR_LIST")

echo "=========================================="
echo "Results: $((TOTAL - NUM_FAILED)) / $TOTAL passed validation"
echo "=========================================="

if [ "$NUM_FAILED" -gt 0 ]; then
    if [ "$ATTEMPT" -ge "$MAX_RETRIES" ]; then
        echo "ERROR: $NUM_FAILED SRRs still failed after $MAX_RETRIES attempts. Giving up."
        echo "Permanently failed SRRs saved to: $FAILED_LIST"
        echo "Manual intervention required for:"
        cat "$FAILED_LIST"
        exit 1
    fi

    NEXT_ATTEMPT=$((ATTEMPT + 1))
    echo "Resubmitting $NUM_FAILED failed SRRs (attempt $NEXT_ATTEMPT)..."
    echo "Failed List: $FAILED_LIST"

    # Submit the rescue download job
    RESCUE_JOB_ID=$(sbatch --parsable --array=1-$NUM_FAILED \
        /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/data_retrieval/download_sra_array.sh \
        "$FAILED_LIST" "$OUT_DIR")

    echo "Submitted rescue download jobs (JobID: $RESCUE_JOB_ID)."

    # Submit next validation sweep AFTER rescue finishes, with incremented attempt counter
    sbatch --dependency=afterany:$RESCUE_JOB_ID \
        /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/data_retrieval/validate_and_rescue.sh \
        "$FAILED_LIST" "$OUT_DIR" "$DATASET_NAME" "$NEXT_ATTEMPT"
else
    echo "=========================================="
    echo "SUCCESS: All $TOTAL SRRs in $DATASET_NAME downloaded and validated!"
    echo "=========================================="
    rm -f "$FAILED_LIST"
fi
