#!/bin/bash
#SBATCH --job-name=dl_gse159911
#SBATCH --output=GSE159911/logs/dl_%A_%a.log
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --array=0-151%50

# Setup
PROJECT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"
GSE_DIR="${PROJECT_DIR}/GSE159911"
METADATA_FILE="${GSE_DIR}/metadata/GSE159911_run_info.csv"
SRA_CACHE="${GSE_DIR}/sra_cache"
FASTQ_DIR="${GSE_DIR}/fastq"

mkdir -p "$SRA_CACHE" "$FASTQ_DIR" "$GSE_DIR/logs"

# Activate Environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl

# Get SRR ID for this task
# Skip header (tail +2), get Nth line (sed), cut 1st column (Run)
# Note: Metadata file is CSV: Run,BioSample
SRR_ID=$(tail -n +2 "$METADATA_FILE" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" | cut -d',' -f1)

if [ -z "$SRR_ID" ]; then
    echo "Error: Could not retrieve SRR ID for index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Processing Task $SLURM_ARRAY_TASK_ID: $SRR_ID"

# 1. Prefetch (Download)
echo "Downloading $SRR_ID..."
prefetch --type sra --output-directory "$SRA_CACHE" "$SRR_ID"

if [ $? -ne 0 ]; then
    echo "Error: Prefetch failed for $SRR_ID"
    exit 1
fi

# 2. Extract (Fasterq-dump)
SRA_FILE="${SRA_CACHE}/${SRR_ID}/${SRR_ID}.sra"
if [ ! -f "$SRA_FILE" ]; then
    # Sometimes prefetch puts it directly in cache root or flat structure
    if [ -f "${SRA_CACHE}/${SRR_ID}.sra" ]; then
        SRA_FILE="${SRA_CACHE}/${SRR_ID}.sra"
    else
        echo "Error: SRA file not found at $SRA_FILE"
        exit 1
    fi
fi

echo "Extracting $SRR_ID..."
fasterq-dump --split-files --threads 4 --outdir "$FASTQ_DIR" --temp "$FASTQ_DIR/temp_${SRR_ID}" "$SRA_FILE"

# 3. Compress (Pigz)
echo "Compressing FASTQs..."
pigz -p 4 "$FASTQ_DIR/${SRR_ID}"*.fastq

# 4. Cleanup
echo "Cleaning up SRA file..."
rm "$SRA_FILE"
rmdir "${SRA_CACHE}/${SRR_ID}" 2>/dev/null

echo "Success: $SRR_ID"
