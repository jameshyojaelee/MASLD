#!/bin/bash
#SBATCH --job-name=download_sra
#SBATCH --output=scripts/data_retrieval/logs/%x_%A_%a.out
#SBATCH --error=scripts/data_retrieval/logs/%x_%A_%a.err
#SBATCH --partition=cpu,io
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00

# Usage: sbatch --array=1-N download_sra_array.sh <SRR_LIST_FILE> <OUTPUT_DIR>

SRR_LIST=$1
OUT_DIR=$2

if [ -z "$SRR_LIST" ] || [ -z "$OUT_DIR" ]; then
    echo "Error: Must provide SRR_LIST and OUT_DIR"
    exit 1
fi

mkdir -p "$OUT_DIR"

# Get the SRR ID for this array task
SRR_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")

if [ -z "$SRR_ID" ]; then
    echo "Error: Could not determine SRR ID for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Processing $SRR_ID..."

# Load SRA Toolkit module (adjust if needed, or assume it's in PATH)
# The lab standard download environment
DL_ENV="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl"
PREFETCH="$DL_ENV/bin/prefetch"
FASTERQ="$DL_ENV/bin/fasterq-dump"
PIGZ="$DL_ENV/bin/pigz"

# Configure SRA toolkit to completely ignore global cache locks across parallel jobs
export VDB_NO_LOCK=1
export NCBI_SETTINGS=/dev/null

# Reroute SRA temporary extraction processes to the high-capacity /scratch volume
export TMPDIR=/scratch
export TEMP=/scratch
export TMP=/scratch

# 1. Prefetch the SRA file
echo "Running prefetch for $SRR_ID..."
$PREFETCH "$SRR_ID" --max-size 100G -O "$OUT_DIR" -f yes

if [ $? -ne 0 ]; then
    echo "Error: prefetch failed for $SRR_ID"
    exit 1
fi

# 2. Convert to FASTQ
echo "Running fasterq-dump for $SRR_ID..."
$FASTERQ --split-files --include-technical --threads 4 -O "$OUT_DIR/$SRR_ID" "$OUT_DIR/$SRR_ID/$SRR_ID.sra" -f

if [ $? -ne 0 ]; then
    echo "Error: fasterq-dump failed for $SRR_ID"
    exit 1
fi

# 3. Compress FASTQ files to save space
echo "Compressing FASTQ files for $SRR_ID..."
$PIGZ -p 4 "$OUT_DIR/$SRR_ID"/*.fastq

# 4. Cleanup raw SRA cache to save quota
echo "Cleaning up .sra for $SRR_ID..."
rm "$OUT_DIR/$SRR_ID/$SRR_ID.sra"

echo "Done $SRR_ID"
