#!/bin/bash
#SBATCH --job-name=liver_multiqc
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=01:00:00
#SBATCH --output=../logs/multiqc/multiqc_%j.out
#SBATCH --error=../logs/multiqc/multiqc_%j.err

# MultiQC report generation for GSE130970 RNA-seq data
# This script generates comprehensive quality control reports
# Requires: Parallel feature counting pipeline (steps 03-04) to be completed

set -euo pipefail

# Source configuration
source ../config.sh

# Load required modules
module load foss/2023b
module load MultiQC/1.22.3-foss-2023b
# module load Python

echo "Starting MultiQC report generation..."
echo "Job ID: ${SLURM_JOB_ID}"
echo "Date: $(date)"

# Create output directory
mkdir -p ${RESULTS_DIR}/multiqc

# Check if parallel feature counting completed successfully
echo "Validating prerequisite steps..."

# Check if count matrix exists (from parallel pipeline)
if [[ ! -f "${RESULTS_DIR}/counts/gene_counts_matrix.txt" ]]; then
    echo "WARNING: Count matrix not found: ${RESULTS_DIR}/counts/gene_counts_matrix.txt"
    echo "MultiQC will generate report with available data (FastQC, STAR, individual counts)"
    echo "Re-run MultiQC after the merge job completes for a complete report"
fi

# Check if merge summary exists
if [[ ! -f "${RESULTS_DIR}/counts/merge_summary.txt" ]]; then
    echo "WARNING: Merge summary not found. Pipeline may not have completed properly."
fi

# Check if individual counts directory exists and has files
if [[ ! -d "${RESULTS_DIR}/counts/individual" ]]; then
    echo "ERROR: Individual counts directory not found: ${RESULTS_DIR}/counts/individual"
    echo "Please run the parallel feature counting pipeline first."
    exit 1
fi

INDIVIDUAL_COUNT_FILES=$(find ${RESULTS_DIR}/counts/individual -name "*_counts.txt" -type f | wc -l)
echo "Found ${INDIVIDUAL_COUNT_FILES} individual count files"

if [[ ${INDIVIDUAL_COUNT_FILES} -eq 0 ]]; then
    echo "WARNING: No individual count files found in ${RESULTS_DIR}/counts/individual"
    echo "MultiQC will generate report with available data (FastQC and STAR only)"
else
    echo "Found ${INDIVIDUAL_COUNT_FILES} individual count files - proceeding with MultiQC"
fi

# Create MultiQC configuration file
cat > ${RESULTS_DIR}/multiqc/multiqc_config.yaml << 'EOF'
title: "GSE130970 NAFLD RNA-seq Analysis - Quality Control Report"
subtitle: "Comprehensive QC report for bulk RNA-seq pipeline"
intro_text: "This report summarizes the quality control metrics for the GSE130970 NAFLD RNA-seq dataset processing pipeline, including FastQC, STAR alignment, and feature counting results."

report_header_info:
    - Contact E-mail: 'jameslee@example.com'
    - Application Type: 'RNA-seq'
    - Project Type: 'NAFLD/NASH Analysis'
    - Sequencing Platform: 'Illumina HiSeq 2500'
    - Read Type: 'Paired-end, 101bp'

extra_fn_clean_exts:
    - '_Aligned.sortedByCoord.out'
    - '_1.fastq'
    - '_2.fastq'
    - '.fastq'

module_order:
    - fastqc:
        name: 'FastQC (Raw Reads)'
        info: 'Quality metrics for raw FASTQ files before any processing.'
    - star:
        name: 'STAR Alignment'
        info: 'Read alignment statistics from STAR aligner.'
    - featurecounts:
        name: 'Feature Counts'
        info: 'Gene quantification results from featureCounts.'
    - samtools:
        name: 'SAMtools Stats'
        info: 'BAM file statistics from SAMtools.'

table_columns_visible:
    FastQC:
        percent_gc: True
        avg_sequence_length: True
        percent_duplicates: True
        total_sequences: True
    STAR:
        uniquely_mapped_percent: True
        multimapped_percent: True
        unmapped_mismatches_percent: True
        unmapped_tooshort_percent: True
    featureCounts:
        percent_assigned: True
        Assigned: True
        Unassigned_MultiMapping: True
        Unassigned_NoFeatures: True

EOF

echo "Searching for analysis files..."

# Check what files are available
echo "Available FastQC files:"
find ${RESULTS_DIR}/fastqc_raw -name "*_fastqc.zip" -type f | wc -l || echo "0"

echo "Available STAR log files:"
find ${RESULTS_DIR}/star_alignment -name "*Log.final.out" -type f | wc -l || echo "0"

echo "Available featureCounts files:"
find ${RESULTS_DIR}/counts/individual -name "*.txt.summary" -type f | wc -l 2>/dev/null || echo "0"

echo "Available SAMtools files:"
find ${RESULTS_DIR}/star_alignment -name "*_alignment_stats.txt" -type f | wc -l || echo "0"

# Run MultiQC
echo "Running MultiQC..."
multiqc \
    --config ${RESULTS_DIR}/multiqc/multiqc_config.yaml \
    --outdir ${RESULTS_DIR}/multiqc \
    --filename GSE130970_multiqc_report \
    --title "GSE130970 NAFLD RNA-seq QC Report" \
    --force \
    --verbose \
    ${RESULTS_DIR}/fastqc_raw \
    ${RESULTS_DIR}/star_alignment \
    ${RESULTS_DIR}/counts/individual

# Check if MultiQC completed successfully
if [[ $? -eq 0 ]]; then
    echo "MultiQC completed successfully"
    
    # Generate additional summary statistics
    echo "Generating additional summary statistics..."
    
    # Create a summary of key metrics
    cat > ${RESULTS_DIR}/multiqc/pipeline_summary.txt << EOF
GSE130970 NAFLD RNA-seq Pipeline Summary
Generated on: $(date)

DATASET INFORMATION:
- Study: GSE130970
- Disease: NAFLD
- Platform: Illumina HiSeq 2500
- Read Type: Paired-end, 101bp
- Total Samples: $(wc -l < ${SAMPLE_LIST})

PROCESSING STEPS COMPLETED:
1. Quality Control (FastQC)
2. Read Alignment (STAR)
3. Gene Quantification (featureCounts - Parallel Processing)
4. Count Matrix Merging
5. Quality Report Generation (MultiQC)

FILES GENERATED:
- MultiQC Report: ${RESULTS_DIR}/multiqc/GSE130970_multiqc_report.html
- Count Matrix: ${RESULTS_DIR}/counts/gene_counts_matrix.txt
- Merge Summary: ${RESULTS_DIR}/counts/merge_summary.txt
- Individual Count Files: ${RESULTS_DIR}/counts/individual/
- Sample Metadata: ${METADATA_FILE}

NEXT STEPS:
- Run differential expression analysis with DESeq2 (script 06_deseq2_analysis.sh)
- Perform functional enrichment analysis
- Generate publication-ready figures

EOF
    
    # Check if we have samples with low mapping rates
    echo "Checking for samples with potential quality issues..."
    if [[ -f "${RESULTS_DIR}/star_alignment/alignment_summary.txt" ]]; then
        grep -E "Uniquely mapped reads %|Sample:" ${RESULTS_DIR}/star_alignment/alignment_summary.txt | \
        paste - - | \
        awk -F'\t|%' '$2 < 70 {print "WARNING: Sample " $1 " has low mapping rate: " $2 "%"}' > ${RESULTS_DIR}/multiqc/quality_warnings.txt
        
        if [[ -s "${RESULTS_DIR}/multiqc/quality_warnings.txt" ]]; then
            echo "Quality warnings found:"
            cat ${RESULTS_DIR}/multiqc/quality_warnings.txt
        else
            echo "All samples pass basic quality thresholds"
        fi
    fi
    
    echo "MultiQC report files:"
    ls -lh ${RESULTS_DIR}/multiqc/GSE130970_multiqc_report.*
    
    echo ""
    echo "=== MULTIQC REPORT READY ==="
    echo "View the report at: ${RESULTS_DIR}/multiqc/GSE130970_multiqc_report.html"
    echo "Pipeline summary: ${RESULTS_DIR}/multiqc/pipeline_summary.txt"
    
else
    echo "ERROR: MultiQC failed"
    exit 1
fi

echo "MultiQC report generation completed at $(date)" 