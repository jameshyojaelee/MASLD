#!/bin/bash
#SBATCH --job-name=validate_proteomics
#SBATCH --output=scripts/data_retrieval/logs/%x_%A.out
#SBATCH --error=scripts/data_retrieval/logs/%x_%A.err
#SBATCH --partition=cpu,io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00

echo "Validating Proteomics Cohorts (GSE276114, PXD052937)..."

FAILED=0

# 1. Validate GSE276114
GSE_FILE="data/Proteomics/GSE276114/GSE276114_raw.count.txt.gz"
echo "Validating $GSE_FILE..."
if [ -f "$GSE_FILE" ]; then
    if ! gzip -t "$GSE_FILE" 2>/dev/null; then
        echo "Corrupted file: $GSE_FILE"
        FAILED=1
    else
        echo "GSE276114 matrix is intact."
    fi
else
    echo "Missing file: $GSE_FILE"
    FAILED=1
fi

# 2. Validate PXD052937
PXD_DIR="data/Proteomics/PXD052937"
# Check if directory exists and has .sne files
if [ ! -d "$PXD_DIR" ] || [ -z "$(find "$PXD_DIR" -type f -name "*.sne")" ]; then
    echo "ERROR: PXD052937 directory is missing or contains no .sne MassSpec matrix files."
    FAILED=1
else
    echo "PXD052937 traces exist. Checking file sizes to ensure completion..."
    # A .sne trace is typically tens of GBs. A file <1MB usually means a broken download link.
    for sne_file in $(find "$PXD_DIR" -type f -name "*.sne"); do
        size=$(stat -c%s "$sne_file")
        if [ "$size" -lt 1000000 ]; then
            echo "Warning: Suspiciously small .sne file detected: $sne_file ($size bytes)"
            FAILED=1
        fi
    done
fi

if [ "$FAILED" -eq 1 ]; then
    echo "------------------------------------------------------"
    echo "Validation FAILED. Corrupt or missing proteomics data found."
    echo "Resubmitting the PXD052937 download script to resume interrupted files via wget -c..."
    echo "------------------------------------------------------"
    
    # Resubmit PXD052937 script
    RESCUE_JOB_ID=$(sbatch --parsable scripts/data_retrieval/download_proteomics_pxd.sh)
    echo "Submitted rescue download job: $RESCUE_JOB_ID"
    
    # Recursively validate again
    sbatch --dependency=afterany:$RESCUE_JOB_ID scripts/data_retrieval/validate_proteomics.sh
else
    echo "------------------------------------------------------"
    echo "SUCCESS: Proteomics files across cohorts GSE276114 and PXD052937 validated successfully."
    echo "------------------------------------------------------"
fi
