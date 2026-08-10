#!/bin/bash
#SBATCH --job-name=fC_merge
#SBATCH --partition=cpu,io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/logs/featurecounts_s2/GSE213621_merge_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/logs/featurecounts_s2/GSE213621_merge_%j.err

echo "ERROR: positional canonical merge is disabled; use bg001_remediation/merge_validate.sbatch." >&2
exit 64

set -euo pipefail

OUT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE213621/counts/featurecounts"
DS="GSE213621"

echo "=== Merging GSE213621 featureCounts chunks ==="
echo "Started: $(date)"

# Use chunk0 as base (has header + annotation columns)
cp "${OUT_DIR}/gene_counts_chunk0.txt" "${OUT_DIR}/gene_counts_merged.txt"

# For chunks 1-3, extract only the count columns (skip first 6 annotation cols) and paste
for i in 1 2 3; do
  CHUNK="${OUT_DIR}/gene_counts_chunk${i}.txt"
  if [[ -f "${CHUNK}" ]]; then
    # Skip header line, cut columns 7+ (sample counts only)
    tail -n +2 "${CHUNK}" | cut -f7- > "${OUT_DIR}/_tmp_chunk${i}_counts.txt"
    # Paste onto merged file (skip header of merged, paste, re-add header)
    head -1 "${OUT_DIR}/gene_counts_merged.txt" > "${OUT_DIR}/_tmp_header.txt"
    head -1 "${CHUNK}" | cut -f7- > "${OUT_DIR}/_tmp_chunk_header.txt"
    paste <(head -1 "${OUT_DIR}/gene_counts_merged.txt") <(cat "${OUT_DIR}/_tmp_chunk_header.txt") > "${OUT_DIR}/_tmp_new_header.txt"
    paste <(tail -n +2 "${OUT_DIR}/gene_counts_merged.txt") <(cat "${OUT_DIR}/_tmp_chunk${i}_counts.txt") > "${OUT_DIR}/_tmp_merged_body.txt"
    cat "${OUT_DIR}/_tmp_new_header.txt" "${OUT_DIR}/_tmp_merged_body.txt" > "${OUT_DIR}/gene_counts_merged.txt"
    rm -f "${OUT_DIR}/_tmp_"*
    echo "  Merged chunk ${i}"
  fi
done

# Replace gene_counts.txt with merged result
mv "${OUT_DIR}/gene_counts_merged.txt" "${OUT_DIR}/gene_counts.txt"
echo "SUCCESS: gene_counts.txt written ($(stat -c%s "${OUT_DIR}/gene_counts.txt") bytes)"
echo "Gene count: $(tail -n +3 "${OUT_DIR}/gene_counts.txt" | wc -l)"

# Cleanup chunks
rm -f "${OUT_DIR}/gene_counts_chunk"*.txt "${OUT_DIR}/gene_counts_chunk"*.txt.summary

echo "Finished: $(date)"
