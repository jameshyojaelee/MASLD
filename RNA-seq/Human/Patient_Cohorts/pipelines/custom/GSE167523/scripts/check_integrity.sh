#!/bin/bash
#SBATCH --job-name=liver_check_integrity
#SBATCH --output=logs/check_%j.out
#SBATCH --error=logs/check_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail

# USAGE: sbatch check_integrity.sh <DIRECTORY_TO_SCAN> <OUTPUT_LOG>

TARGET_DIR=$1
LOG_FILE=$2

echo "Scanning $TARGET_DIR..."
rm -f "$LOG_FILE"

find "$TARGET_DIR" -name "*.fastq.gz" | while read -r FILE; do
    if ! gzip -t "$FILE"; then
        echo "CORRUPT: $FILE"
        echo "$FILE" >> "$LOG_FILE"
    else
        echo "OK: $FILE"
    fi
done

echo "Scan complete. Corrupt files listed in $LOG_FILE"
