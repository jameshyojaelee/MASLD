#!/bin/bash
#SBATCH --job-name=liver_star
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60GB
#SBATCH --time=06:00:00
#SBATCH --array=1-79%4
#SBATCH --output=../logs/star/star_%A_%a.out
#SBATCH --error=../logs/star/star_%A_%a.err

# STAR alignment for GSE130970 RNA-seq data
# This script aligns paired-end reads to the reference genome using STAR

set -euo pipefail

# Source configuration
source ../config.sh

# Load required modules
module load STAR/2.7.11b-GCC-13.2.0
module load SAMtools/1.21

echo "Starting STAR alignment..."
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

# Create output directory
SAMPLE_OUT_DIR="${RESULTS_DIR}/star_alignment/${SAMPLE}"
mkdir -p ${SAMPLE_OUT_DIR}

# Check if STAR index exists
if [[ ! -d "${STAR_INDEX}" ]]; then
    echo "ERROR: STAR index not found at ${STAR_INDEX}"
    exit 1
fi

echo "Input files:"
echo "R1: ${R1} ($(ls -lh ${R1} | awk '{print $5}'))"
echo "R2: ${R2} ($(ls -lh ${R2} | awk '{print $5}'))"
echo "Output directory: ${SAMPLE_OUT_DIR}"

# Run STAR alignment
echo "Running STAR alignment for ${SAMPLE}..."
STAR \
    --runThreadN ${STAR_THREADS} \
    --genomeDir ${STAR_INDEX} \
    --readFilesIn ${R1} ${R2} \
    --outFileNamePrefix ${SAMPLE_OUT_DIR}/${SAMPLE}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --quantMode GeneCounts \
    --sjdbOverhang 100 \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonical \
    --chimSegmentMin 20 \
    --chimJunctionOverhangMin 20 \
    --chimOutType Junctions SeparateSAMold WithinBAM \
    --chimMainSegmentMultNmax 1

# Check if STAR completed successfully
if [[ $? -eq 0 ]]; then
    echo "STAR alignment completed successfully for ${SAMPLE}"
    
    # Index the BAM file
    echo "Indexing BAM file..."
    samtools index ${SAMPLE_OUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam
    
    # Generate alignment statistics
    echo "Generating alignment statistics..."
    samtools flagstat ${SAMPLE_OUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam > ${SAMPLE_OUT_DIR}/${SAMPLE}_alignment_stats.txt
    samtools idxstats ${SAMPLE_OUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam > ${SAMPLE_OUT_DIR}/${SAMPLE}_idxstats.txt
    
    # Extract key statistics from STAR log
    LOG_FILE="${SAMPLE_OUT_DIR}/${SAMPLE}_Log.final.out"
    if [[ -f "${LOG_FILE}" ]]; then
        echo "Sample: ${SAMPLE}" >> ${RESULTS_DIR}/star_alignment/alignment_summary.txt
        grep "Number of input reads" ${LOG_FILE} >> ${RESULTS_DIR}/star_alignment/alignment_summary.txt
        grep "Uniquely mapped reads number" ${LOG_FILE} >> ${RESULTS_DIR}/star_alignment/alignment_summary.txt
        grep "Uniquely mapped reads %" ${LOG_FILE} >> ${RESULTS_DIR}/star_alignment/alignment_summary.txt
        grep "Number of reads mapped to multiple loci" ${LOG_FILE} >> ${RESULTS_DIR}/star_alignment/alignment_summary.txt
        grep "% of reads mapped to multiple loci" ${LOG_FILE} >> ${RESULTS_DIR}/star_alignment/alignment_summary.txt
        grep "Number of reads unmapped" ${LOG_FILE} >> ${RESULTS_DIR}/star_alignment/alignment_summary.txt
        grep "% of reads unmapped" ${LOG_FILE} >> ${RESULTS_DIR}/star_alignment/alignment_summary.txt
        echo "Processed at: $(date)" >> ${RESULTS_DIR}/star_alignment/alignment_summary.txt
        echo "---" >> ${RESULTS_DIR}/star_alignment/alignment_summary.txt
    fi
    
    # Cleanup temporary files to save space
    echo "Cleaning up temporary files..."
    rm -f ${SAMPLE_OUT_DIR}/${SAMPLE}_Aligned.out.sam 2>/dev/null || true
    rm -rf ${SAMPLE_OUT_DIR}/${SAMPLE}__STARtmp 2>/dev/null || true
    
    echo "File sizes after alignment:"
    ls -lh ${SAMPLE_OUT_DIR}/${SAMPLE}_*.bam ${SAMPLE_OUT_DIR}/${SAMPLE}_*.txt 2>/dev/null || true
    
else
    echo "ERROR: STAR alignment failed for ${SAMPLE}"
    exit 1
fi

echo "STAR alignment completed for ${SAMPLE} at $(date)" 