#!/bin/bash
#SBATCH --job-name=liver_repair_corruption
#SBATCH --output=logs/repair_%j.out
#SBATCH --error=logs/repair_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

# 1. Setup Environment
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

# Function: Download and Compress
# Usage: download_sample <ACCESSION> <DEST_DIR>
download_sample() {
    local ACC=$1
    local DEST=$2
    
    echo "[start] Repairing $ACC -> $DEST"
    mkdir -p "$DEST"
    
    fasterq-dump --outdir "$DEST" --split-files --threads 8 --force --progress "$ACC"
    gzip -f "$DEST"/${ACC}*.fastq
    echo "[done] $ACC repaired."
}

# 2. Repair GSE167523 Sample (SRR13797146)
download_sample "SRR13797146" "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE167523/fastq"

# 3. Repair GSE126848 Sample (SRR8601578)
download_sample "SRR8601578" "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE126848/fastq"

echo "All repairs complete."
