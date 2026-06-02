#!/bin/bash
#SBATCH --job-name=liver_merge_counts
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=01:00:00
#SBATCH --output=../logs/featurecounts/merge_counts_%j.out
#SBATCH --error=../logs/featurecounts/merge_counts_%j.err

# Merge individual feature count files into a single count matrix
# Run this after all parallel feature counting jobs complete

set -euo pipefail

# Source configuration
source ../config.sh

echo "Starting count matrix merge..."
echo "Job ID: ${SLURM_JOB_ID}"
echo "Date: $(date)"

# Check that counts directory exists
if [[ ! -d "${RESULTS_DIR}/counts/individual" ]]; then
    echo "ERROR: Individual counts directory not found: ${RESULTS_DIR}/counts/individual"
    exit 1
fi

# Check completion status of all samples
echo "Checking completion status of all samples..."
TOTAL_SAMPLES=0
COMPLETED_SAMPLES=0
FAILED_SAMPLES=0

while read sample; do
    ((TOTAL_SAMPLES++))
    STATUS_FILE="${RESULTS_DIR}/counts/individual/${sample}_status.txt"
    COUNT_FILE="${RESULTS_DIR}/counts/individual/${sample}_counts.txt"
    
    if [[ -f "${STATUS_FILE}" && "$(cat ${STATUS_FILE})" == "COMPLETED" && -f "${COUNT_FILE}" ]]; then
        ((COMPLETED_SAMPLES++))
    else
        echo "WARNING: Sample ${sample} not completed or missing files"
        ((FAILED_SAMPLES++))
    fi
done < ${SAMPLE_LIST}

echo "Sample status: ${COMPLETED_SAMPLES} completed, ${FAILED_SAMPLES} failed/missing out of ${TOTAL_SAMPLES} total"

if [[ ${COMPLETED_SAMPLES} -eq 0 ]]; then
    echo "ERROR: No samples were processed successfully"
    exit 1
fi

# Get the first successfully completed file to extract gene IDs
echo "Finding reference file for gene IDs..."
FIRST_FILE=""
while read sample; do
    STATUS_FILE="${RESULTS_DIR}/counts/individual/${sample}_status.txt"
    COUNT_FILE="${RESULTS_DIR}/counts/individual/${sample}_counts.txt"
    
    if [[ -f "${STATUS_FILE}" && "$(cat ${STATUS_FILE})" == "COMPLETED" && -f "${COUNT_FILE}" ]]; then
        FIRST_FILE="${COUNT_FILE}"
        echo "Using ${sample} as reference file"
        break
    fi
done < ${SAMPLE_LIST}

if [[ -z "${FIRST_FILE}" ]]; then
    echo "ERROR: No completed count files found"
    exit 1
fi

# Extract gene IDs (first column, skip header)
echo "Extracting gene IDs..."
tail -n +3 "${FIRST_FILE}" | cut -f1 > "${RESULTS_DIR}/counts/gene_ids.txt"
GENE_COUNT=$(wc -l < "${RESULTS_DIR}/counts/gene_ids.txt")
echo "Found ${GENE_COUNT} genes"

# Create the header for the final matrix
echo "Building count matrix header..."
echo -n "GeneID" > "${RESULTS_DIR}/counts/gene_counts_matrix.txt"

# Add sample names to header (only for completed samples)
INCLUDED_SAMPLES=0
while read sample; do
    STATUS_FILE="${RESULTS_DIR}/counts/individual/${sample}_status.txt"
    COUNT_FILE="${RESULTS_DIR}/counts/individual/${sample}_counts.txt"
    
    if [[ -f "${STATUS_FILE}" && "$(cat ${STATUS_FILE})" == "COMPLETED" && -f "${COUNT_FILE}" ]]; then
        echo -ne "\t${sample}" >> "${RESULTS_DIR}/counts/gene_counts_matrix.txt"
        ((INCLUDED_SAMPLES++))
    fi
done < ${SAMPLE_LIST}
echo "" >> "${RESULTS_DIR}/counts/gene_counts_matrix.txt"

echo "Matrix will include ${INCLUDED_SAMPLES} samples"

# Extract counts for each completed sample and create temporary files
echo "Extracting count data..."
TEMP_FILES=()
TEMP_COUNTER=0

while read sample; do
    STATUS_FILE="${RESULTS_DIR}/counts/individual/${sample}_status.txt"
    COUNT_FILE="${RESULTS_DIR}/counts/individual/${sample}_counts.txt"
    
    if [[ -f "${STATUS_FILE}" && "$(cat ${STATUS_FILE})" == "COMPLETED" && -f "${COUNT_FILE}" ]]; then
        # Extract just the count column (column 7) and save to temp file
        TEMP_FILE="${RESULTS_DIR}/counts/individual/${sample}_counts_only.txt"
        tail -n +3 "${COUNT_FILE}" | cut -f7 > "${TEMP_FILE}"
        TEMP_FILES+=("${TEMP_FILE}")
        ((TEMP_COUNTER++))
        
        if [[ $((TEMP_COUNTER % 20)) -eq 0 ]]; then
            echo "Processed ${TEMP_COUNTER}/${INCLUDED_SAMPLES} samples"
        fi
    fi
done < ${SAMPLE_LIST}

echo "Merging count data..."
# Paste gene IDs with all count columns
paste "${RESULTS_DIR}/counts/gene_ids.txt" "${TEMP_FILES[@]}" >> "${RESULTS_DIR}/counts/gene_counts_matrix.txt"

# Verify the final matrix
MATRIX_ROWS=$(wc -l < "${RESULTS_DIR}/counts/gene_counts_matrix.txt")
MATRIX_COLS=$(head -n1 "${RESULTS_DIR}/counts/gene_counts_matrix.txt" | wc -w)

echo "Count matrix created: ${RESULTS_DIR}/counts/gene_counts_matrix.txt"
echo "Matrix dimensions: ${MATRIX_ROWS} rows (including header), ${MATRIX_COLS} columns"
echo "Expected: $((GENE_COUNT + 1)) rows, $((INCLUDED_SAMPLES + 1)) columns"

# Create a summary file
SUMMARY_FILE="${RESULTS_DIR}/counts/merge_summary.txt"
echo "Count Matrix Merge Summary" > "${SUMMARY_FILE}"
echo "Date: $(date)" >> "${SUMMARY_FILE}"
echo "Total samples: ${TOTAL_SAMPLES}" >> "${SUMMARY_FILE}"
echo "Completed samples: ${COMPLETED_SAMPLES}" >> "${SUMMARY_FILE}"
echo "Failed samples: ${FAILED_SAMPLES}" >> "${SUMMARY_FILE}"
echo "Samples included in matrix: ${INCLUDED_SAMPLES}" >> "${SUMMARY_FILE}"
echo "Genes: ${GENE_COUNT}" >> "${SUMMARY_FILE}"
echo "Matrix file: gene_counts_matrix.txt" >> "${SUMMARY_FILE}"

# List failed samples if any
if [[ ${FAILED_SAMPLES} -gt 0 ]]; then
    echo "" >> "${SUMMARY_FILE}"
    echo "Failed/missing samples:" >> "${SUMMARY_FILE}"
    while read sample; do
        STATUS_FILE="${RESULTS_DIR}/counts/individual/${sample}_status.txt"
        COUNT_FILE="${RESULTS_DIR}/counts/individual/${sample}_counts.txt"
        
        if [[ ! -f "${STATUS_FILE}" || "$(cat ${STATUS_FILE} 2>/dev/null)" != "COMPLETED" || ! -f "${COUNT_FILE}" ]]; then
            echo "  ${sample}" >> "${SUMMARY_FILE}"
        fi
    done < ${SAMPLE_LIST}
fi

# Clean up temporary files
echo "Cleaning up temporary files..."
rm -f "${RESULTS_DIR}/counts/gene_ids.txt"
rm -f "${TEMP_FILES[@]}"

echo "Count matrix merge completed successfully at $(date)"
echo "Summary written to: ${SUMMARY_FILE}" 