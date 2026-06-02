#!/bin/bash
#SBATCH --job-name=liver_integration_6cohort
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/6cohort_integration_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/6cohort_integration_%j.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --time=12:00:00
# NOTE: dependency is managed by merge_featurecounts_and_launch.sh — do NOT add --dependency here

# =============================================================================
# 6-Cohort Integration Pipeline with GSE213621 Safety Checks
# =============================================================================
# This script:
#   1. Validates GSE213621 gene_counts.txt integrity and completeness
#   2. Runs the FULL integration pipeline (scripts 00 → 12) on all 6 cohorts
#
# Submit with:
#   sbatch run_6cohort_integration_slurm.sh
#
# The --dependency=afterok:13880407 ensures this job only starts after the
# GSE213621 Snakemake pipeline (submit_pipeline) completes successfully.
# =============================================================================

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "============================================================"
echo "  6-COHORT MASLD INTEGRATION PIPELINE"
echo "  Started: $(date)"
echo "  SLURM_JOB_ID: $SLURM_JOB_ID"
echo "  Partition: bigmem | RAM: 400G | CPUs: $SLURM_CPUS_PER_TASK"
echo "============================================================"
echo ""

# =============================================================================
# PHASE 0: SAFETY CHECKS — Validate GSE213621 gene_counts.txt
# =============================================================================

COUNTS_FILE="results/GSE213621/counts/featurecounts/gene_counts.txt"
COUNTS_SUMMARY="${COUNTS_FILE}.summary"
EXPECTED_SAMPLES=367    # 367 data samples from SraRunTable (368 rows - 1 header)
EXPECTED_COLS=$((EXPECTED_SAMPLES + 6))  # 6 featureCounts annotation cols + samples
MIN_GENES=30000         # GENCODE v49 has ~60k+ entries; minimum sanity threshold
MIN_FILE_SIZE_MB=100    # Batch-merged featureCounts produces ~175MB for 367 samples

FAIL=0

echo "=== PHASE 0: Safety Checks ==="
echo ""

# Check 1: File existence
echo -n "[CHECK 1] gene_counts.txt exists ... "
if [ ! -f "$COUNTS_FILE" ]; then
    echo "FAIL — file not found: $COUNTS_FILE"
    echo "The Snakemake pipeline may not have completed successfully."
    exit 1
fi
echo "PASS"

# Check 2: File is not empty and meets minimum size
FILE_SIZE_BYTES=$(stat -c%s "$COUNTS_FILE")
FILE_SIZE_MB=$((FILE_SIZE_BYTES / 1024 / 1024))
echo -n "[CHECK 2] File size >= ${MIN_FILE_SIZE_MB}MB ... "
if [ "$FILE_SIZE_MB" -lt "$MIN_FILE_SIZE_MB" ]; then
    echo "FAIL — file is only ${FILE_SIZE_MB}MB (expected >= ${MIN_FILE_SIZE_MB}MB)"
    echo "The file may be truncated or incomplete."
    FAIL=1
else
    echo "PASS (${FILE_SIZE_MB}MB)"
fi

# Check 3: File is not being actively written (no change in 30 seconds)
echo -n "[CHECK 3] File not actively being written ... "
SIZE_A=$(stat -c%s "$COUNTS_FILE")
sleep 30
SIZE_B=$(stat -c%s "$COUNTS_FILE")
if [ "$SIZE_A" != "$SIZE_B" ]; then
    echo "FAIL — file size changed from $SIZE_A to $SIZE_B in 30 seconds"
    echo "featureCounts may still be writing. Aborting."
    exit 1
fi
echo "PASS (stable)"

# Check 4: Correct number of columns (samples)
ACTUAL_COLS=$(head -2 "$COUNTS_FILE" | tail -1 | awk -F'\t' '{print NF}')
echo -n "[CHECK 4] Column count = $EXPECTED_COLS (6 anno + $EXPECTED_SAMPLES samples) ... "
if [ "$ACTUAL_COLS" -ne "$EXPECTED_COLS" ]; then
    echo "FAIL — found $ACTUAL_COLS columns (expected $EXPECTED_COLS)"
    echo "Sample count mismatch. Check if all BAMs were included in featureCounts."
    FAIL=1
else
    echo "PASS ($ACTUAL_COLS columns)"
fi

# Check 5: Sufficient gene rows
# Subtract 2 for the comment header line and the column header line
TOTAL_LINES=$(wc -l < "$COUNTS_FILE")
GENE_ROWS=$((TOTAL_LINES - 2))
echo -n "[CHECK 5] Gene count >= $MIN_GENES ... "
if [ "$GENE_ROWS" -lt "$MIN_GENES" ]; then
    echo "FAIL — only $GENE_ROWS genes (expected >= $MIN_GENES)"
    FAIL=1
else
    echo "PASS ($GENE_ROWS genes)"
fi

# Check 6: Header format is valid featureCounts
echo -n "[CHECK 6] Valid featureCounts header ... "
FIRST_LINE=$(head -1 "$COUNTS_FILE")
if [[ "$FIRST_LINE" != "# Program:featureCounts"* ]]; then
    echo "FAIL — first line does not start with '# Program:featureCounts'"
    FAIL=1
else
    echo "PASS"
fi

# Check 7: First data column is Geneid
echo -n "[CHECK 7] First column header is 'Geneid' ... "
FIRST_COL=$(head -2 "$COUNTS_FILE" | tail -1 | cut -f1)
if [ "$FIRST_COL" != "Geneid" ]; then
    echo "FAIL — first column is '$FIRST_COL', not 'Geneid'"
    FAIL=1
else
    echo "PASS"
fi

# Check 8: Last row is not truncated (has expected number of fields)
echo -n "[CHECK 8] Last row integrity ... "
LAST_ROW_COLS=$(tail -1 "$COUNTS_FILE" | awk -F'\t' '{print NF}')
if [ "$LAST_ROW_COLS" -ne "$ACTUAL_COLS" ]; then
    echo "FAIL — last row has $LAST_ROW_COLS columns (header has $ACTUAL_COLS)"
    echo "File may be truncated."
    FAIL=1
else
    echo "PASS"
fi

# Check 9: No zero-only sample columns (all zeroes = failed alignment)
echo -n "[CHECK 9] No dead sample columns (all zeros) ... "
ZERO_COLS=$(awk -F'\t' 'NR>2 {for(i=7;i<=NF;i++) sum[i]+=$i} END {c=0; for(i in sum) if(sum[i]==0) c++; print c}' "$COUNTS_FILE")
if [ "$ZERO_COLS" -gt 0 ]; then
    echo "WARNING — $ZERO_COLS sample column(s) have zero total counts"
    echo "These samples may have failed alignment. Proceeding but flagging."
else
    echo "PASS (no dead columns)"
fi

# Check 10: Summary file exists and is consistent
echo -n "[CHECK 10] featureCounts summary exists ... "
if [ ! -f "$COUNTS_SUMMARY" ]; then
    echo "WARNING — $COUNTS_SUMMARY not found (non-critical)"
else
    SUMMARY_SAMPLES=$(head -1 "$COUNTS_SUMMARY" | awk -F'\t' '{print NF-1}')
    if [ "$SUMMARY_SAMPLES" -ne "$EXPECTED_SAMPLES" ]; then
        echo "WARNING — summary has $SUMMARY_SAMPLES samples, expected $EXPECTED_SAMPLES"
    else
        echo "PASS ($SUMMARY_SAMPLES samples in summary)"
    fi
fi

echo ""

# --- Verdict ---
if [ "$FAIL" -gt 0 ]; then
    echo "╔══════════════════════════════════════════════════╗"
    echo "║  SAFETY CHECK FAILED — ABORTING INTEGRATION     ║"
    echo "║  Review the errors above before proceeding.     ║"
    echo "╚══════════════════════════════════════════════════╝"
    exit 1
fi

echo "╔══════════════════════════════════════════════════╗"
echo "║  ALL SAFETY CHECKS PASSED                       ║"
echo "║  Proceeding with 6-cohort integration pipeline  ║"
echo "╚══════════════════════════════════════════════════╝"
echo ""

# =============================================================================
# PHASE 1: Re-run metadata harmonization (adds GSE213621 to unified_metadata)
# =============================================================================
echo "=== 00: Harmonize Metadata ==="
Rscript analysis/integration/scripts/00_harmonize_metadata.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 00_harmonize_metadata.R"; exit 1; fi

# =============================================================================
# PHASE 2: QC + Per-study DE (re-runs with GSE213621 included)
# =============================================================================
echo ""
echo "=== 01: Sample QC ==="
Rscript analysis/integration/scripts/01_sample_qc.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 01_sample_qc.R"; exit 1; fi

echo ""
echo "=== 02: Per-Study DE ==="
Rscript analysis/integration/scripts/02_per_study_de.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 02_per_study_de.R"; exit 1; fi

# =============================================================================
# PHASE 3: Integration + mega-analysis (the core)
# =============================================================================
echo ""
echo "=== 03: Integrate Counts ==="
Rscript analysis/integration/scripts/03_integrate_counts.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 03_integrate_counts.R"; exit 1; fi

echo ""
echo "=== 04: Variance Partition ==="
Rscript analysis/integration/scripts/04_variance_partition.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 04_variance_partition.R"; exit 1; fi

echo ""
echo "=== 05: Dream Mega-Analysis ($SLURM_CPUS_PER_TASK CPUs) ==="
Rscript analysis/integration/scripts/05_dream_mega_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 05_dream_mega_analysis.R"; exit 1; fi

# =============================================================================
# PHASE 4: Downstream analysis
# =============================================================================
echo ""
echo "=== 07: Consensus DEGs ==="
Rscript analysis/integration/scripts/07_consensus_degs.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 07_consensus_degs.R"; exit 1; fi

echo ""
echo "=== 08: Pathway Analysis ==="
Rscript analysis/integration/scripts/08_pathway_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 08_pathway_analysis.R"; exit 1; fi

echo ""
echo "=== 09: Batch Correction + UMAP ==="
Rscript analysis/integration/scripts/09_batch_correction_umap.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 09_batch_correction_umap.R"; exit 1; fi

echo ""
echo "=== 10: Volcano Plots ==="
Rscript analysis/integration/scripts/10_volcano_plots.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 10_volcano_plots.R"; exit 1; fi

echo ""
echo "=== 12: Library Intersection ==="
Rscript analysis/integration/scripts/12_library_intersection.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 12_library_intersection.R"; exit 1; fi

echo ""
echo "============================================================"
echo "  6-COHORT INTEGRATION COMPLETE"
echo "  Finished: $(date)"
echo "  All scripts (00-12) passed successfully."
echo "  Total samples: ~1,000+ (5 existing + 367 GSE213621)"
echo "============================================================"
