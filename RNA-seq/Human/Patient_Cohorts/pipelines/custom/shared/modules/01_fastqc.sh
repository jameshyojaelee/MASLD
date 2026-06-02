#!/bin/bash
#SBATCH --job-name=liver_fastqc
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --array=1-79%10
#SBATCH --output=../logs/fastqc/fastqc_%A_%a.out
#SBATCH --error=../logs/fastqc/fastqc_%A_%a.err

# FastQC Quality Control for GSE130970 RNA-seq data
# This script runs FastQC on raw FASTQ files using job arrays

set -euo pipefail

# Source configuration
source ../config.sh

# Load required modules
module load fastqc/0.12.1

echo "Starting FastQC analysis..."
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"

# Get sample name from array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})
echo "Processing sample: ${SAMPLE}"

# Define input files
R1="${FASTQ_DIR}/${SAMPLE}_1.fastq"
R2="${FASTQ_DIR}/${SAMPLE}_2.fastq"

# Check if input files exist
if [[ ! -f "$R1" ]] || [[ ! -f "$R2" ]]; then
    echo "ERROR: FASTQ files not found for sample ${SAMPLE}"
    echo "Expected: $R1"
    echo "Expected: $R2"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p ${RESULTS_DIR}/fastqc_raw

# Run FastQC
echo "Running FastQC on ${SAMPLE}..."
fastqc \
    --outdir ${RESULTS_DIR}/fastqc_raw \
    --threads 4 \
    --format fastq \
    --quiet \
    ${R1} ${R2}

# Check if FastQC completed successfully
if [[ $? -eq 0 ]]; then
    echo "FastQC completed successfully for ${SAMPLE}"
    
    # Generate summary statistics
    echo "Sample: ${SAMPLE}" >> ${RESULTS_DIR}/fastqc_raw/fastqc_summary.txt
    echo "R1 size: $(ls -lh ${R1} | awk '{print $5}')" >> ${RESULTS_DIR}/fastqc_raw/fastqc_summary.txt
    echo "R2 size: $(ls -lh ${R2} | awk '{print $5}')" >> ${RESULTS_DIR}/fastqc_raw/fastqc_summary.txt
    echo "Processed at: $(date)" >> ${RESULTS_DIR}/fastqc_raw/fastqc_summary.txt
    echo "---" >> ${RESULTS_DIR}/fastqc_raw/fastqc_summary.txt
else
    echo "ERROR: FastQC failed for ${SAMPLE}"
    exit 1
fi

echo "FastQC analysis completed for ${SAMPLE} at $(date)" 