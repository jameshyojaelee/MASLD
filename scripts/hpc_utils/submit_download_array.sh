#!/bin/bash
#SBATCH --job-name=sra_dl_array
#SBATCH --partition=io
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=4:00:00
#SBATCH --array=1-1000%20
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq/logs/dl_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq/logs/dl_%A_%a.err

# Directory Setup
BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"
LOG_DIR="$BASE_DIR/logs"
job_list="$BASE_DIR/master_download_list.txt"

mkdir -p "$LOG_DIR"

# Environment
source /gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash 2>/dev/null
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl

# Get job details from line N
if [ ! -f "$job_list" ]; then
    echo "Error: Job list $job_list not found!"
    exit 1
fi

# Extract params: GSE_ID, SRR_ID
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$job_list")
if [ -z "$line" ]; then
    echo "No job for index $SLURM_ARRAY_TASK_ID"
    exit 0
fi

GSE_ID=$(echo "$line" | awk '{print $1}')
SRR_ID=$(echo "$line" | awk '{print $2}')

TARGET_DIR="$BASE_DIR/$GSE_ID/fastq"
TEMP_DIR="$BASE_DIR/$GSE_ID/sra_cache"

mkdir -p "$TARGET_DIR" "$TEMP_DIR"

echo "Processing $GSE_ID : $SRR_ID on $(hostname)"

# Check if done
if [ -f "$TARGET_DIR/${SRR_ID}_1.fastq.gz" ] && [ -f "$TARGET_DIR/${SRR_ID}_2.fastq.gz" ]; then
    echo "Files already exist. Skipping..."
    exit 0
fi

# Check if uncompressed files exist (Resume)
if compgen -G "$TARGET_DIR/${SRR_ID}*.fastq" > /dev/null; then
    echo "Found uncompressed FASTQ files. Skipping download/extraction..."
else
    # 1. Prefetch (Download .sra)
    echo "Prefetching..."
    prefetch "$SRR_ID" --output-directory "$TEMP_DIR" --max-size 100G
    if [ $? -ne 0 ]; then
        echo "Prefetch failed."
        exit 1
    fi

    # 2. Extract (fasterq-dump)
    echo "Extracting..."
    fasterq-dump "$TEMP_DIR/$SRR_ID/$SRR_ID.sra" --split-files --outdir "$TARGET_DIR" --threads 4 --include-technical
    if [ $? -ne 0 ]; then
        echo "Extraction failed."
        rm -rf "$TEMP_DIR/$SRR_ID"
        exit 1
    fi
    
    # Cleanup SRA immediately to save space
    rm -rf "$TEMP_DIR/$SRR_ID"
fi

# 3. Compress (pigz)
# Only run if uncompressed files exist
if compgen -G "$TARGET_DIR/${SRR_ID}*.fastq" > /dev/null; then
    echo "Compressing..."
    pigz -p 4 -f "$TARGET_DIR/${SRR_ID}"*.fastq
    if [ $? -ne 0 ]; then
        echo "Compression failed."
        exit 1
    fi
else
    echo "No FASTQ files found to compress!"
    exit 1
fi

echo "Success!"
