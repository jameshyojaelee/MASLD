#!/bin/bash
#SBATCH --job-name=liver_fc_merge_launch
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/slurm_logs/fc_merge_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/slurm_logs/fc_merge_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=2:00:00

echo "ERROR: merge-and-auto-launch is disabled; use bg001_remediation/merge_validate.sbatch and the separate analysis approval gate." >&2
exit 64

# =============================================================================
# Merge parallel featureCounts batches → gene_counts.txt
# Then submit the remaining Snakemake steps (multiqc + deseq2) and the
# 6-cohort integration pipeline.
# =============================================================================

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

export OUTDIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE213621/counts/featurecounts"
export BATCH_DIR="$OUTDIR/batches"
export FINAL_OUT="$OUTDIR/gene_counts.txt"
export FINAL_SUMMARY="$OUTDIR/gene_counts.txt.summary"
export EXPECTED_SAMPLES=367

echo "=== Merging featureCounts Batches ==="
echo "Started: $(date)"
echo ""

# --- Validate all batch files exist ---
FAIL=0
for i in $(seq 0 19); do
    BATCH_FILE="$BATCH_DIR/batch_${i}_counts.txt"
    if [ ! -f "$BATCH_FILE" ]; then
        echo "MISSING: $BATCH_FILE"
        FAIL=1
    fi
done

if [ "$FAIL" -gt 0 ]; then
    echo "ERROR: Some batch files are missing. Aborting merge."
    exit 1
fi
echo "All 20 batch files present."

# --- Merge using Python (handles column alignment properly) ---
python3 << 'PYEOF'
import os, sys

batch_dir = os.environ["BATCH_DIR"]
final_out = os.environ["FINAL_OUT"]
final_summary = os.environ["FINAL_SUMMARY"]
expected_samples = int(os.environ["EXPECTED_SAMPLES"])

# featureCounts output format:
# Line 1: # Program:featureCounts ... (comment header with the command)
# Line 2: Geneid\tChr\tStart\tEnd\tStrand\tLength\t<BAM_path_1>\t<BAM_path_2>\t...
# Line 3+: data rows

anno_cols = 6  # Geneid, Chr, Start, End, Strand, Length

# Read the first batch to get the comment header & annotation columns
first_batch = os.path.join(batch_dir, "batch_0_counts.txt")
with open(first_batch) as f:
    comment_header = f.readline()  # Line 1: # Program:featureCounts ...
    col_header = f.readline().rstrip("\n").split("\t")

# Annotation columns from the first batch
anno_headers = col_header[:anno_cols]

# Collect all sample columns across batches
all_sample_headers = []
# Store data: gene_id -> [anno_fields..., count1, count2, ...]
gene_data = {}  # Geneid -> annotation fields
gene_counts = {}  # Geneid -> list of counts across all batches

for batch_id in range(20):
    batch_file = os.path.join(batch_dir, f"batch_{batch_id}_counts.txt")
    print(f"Reading {batch_file}...", flush=True)
    
    with open(batch_file) as f:
        f.readline()  # skip comment
        header = f.readline().rstrip("\n").split("\t")
        sample_headers = header[anno_cols:]
        all_sample_headers.extend(sample_headers)
        
        for line in f:
            fields = line.rstrip("\n").split("\t")
            gene_id = fields[0]
            anno = fields[:anno_cols]
            counts = fields[anno_cols:]
            
            if gene_id not in gene_data:
                gene_data[gene_id] = anno
                gene_counts[gene_id] = []
            
            gene_counts[gene_id].extend(counts)

total_samples = len(all_sample_headers)
print(f"\nTotal samples merged: {total_samples}")

if total_samples != expected_samples:
    print(f"WARNING: Expected {expected_samples} samples, got {total_samples}")

# Write merged output
print(f"Writing {final_out}...", flush=True)
with open(final_out, "w") as f:
    f.write(comment_header)  # Keep original comment header
    f.write("\t".join(anno_headers + all_sample_headers) + "\n")
    
    # Write in the same order as first batch (preserves gene ordering)
    first_batch_genes = []
    with open(os.path.join(batch_dir, "batch_0_counts.txt")) as fb:
        fb.readline()  # comment
        fb.readline()  # header
        for line in fb:
            first_batch_genes.append(line.split("\t")[0])
    
    for gene_id in first_batch_genes:
        anno = gene_data[gene_id]
        counts = gene_counts[gene_id]
        f.write("\t".join(anno + counts) + "\n")

# Verify
num_genes = len(first_batch_genes)
print(f"Written: {num_genes} genes × {total_samples} samples")
print(f"Total columns: {anno_cols + total_samples}")

# Merge summary files
print(f"\nMerging summary files...", flush=True)
summary_data = {}
summary_headers = []
for batch_id in range(20):
    summary_file = os.path.join(batch_dir, f"batch_{batch_id}_counts.txt.summary")
    if not os.path.exists(summary_file):
        print(f"WARNING: {summary_file} not found, skipping summary merge")
        continue
    with open(summary_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        sample_cols = header[1:]
        summary_headers.extend(sample_cols)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            status = fields[0]
            values = fields[1:]
            if status not in summary_data:
                summary_data[status] = []
            summary_data[status].extend(values)

if summary_headers:
    with open(final_summary, "w") as f:
        f.write("Status\t" + "\t".join(summary_headers) + "\n")
        for status in sorted(summary_data.keys()):
            f.write(status + "\t" + "\t".join(summary_data[status]) + "\n")
    print(f"Summary written: {final_summary}")

print("\n=== Merge complete ===")
PYEOF

if [ $? -ne 0 ]; then
    echo "MERGE FAILED"
    exit 1
fi

echo ""
echo "=== Verifying merged gene_counts.txt ==="
echo "File size: $(ls -lh "$FINAL_OUT" | awk '{print $5}')"
echo "Columns: $(head -2 "$FINAL_OUT" | tail -1 | awk -F'\t' '{print NF}')"
echo "Rows: $(wc -l < "$FINAL_OUT")"
echo ""

# --- Now run multiqc (the other pending Snakemake rule) ---
echo "=== Running MultiQC ==="
RESULTS_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE213621"
MULTIQC_DIR="$RESULTS_ROOT/qc/multiqc"
mkdir -p "$MULTIQC_DIR"
multiqc "$RESULTS_ROOT" -o "$MULTIQC_DIR" 2>&1 || echo "MultiQC warning (non-fatal)"

echo ""
echo "=== Submitting 6-Cohort Integration Pipeline ==="

# Submit the integration job (no dependency needed — we're already past the gate)
INTEGRATION_JOB=$(sbatch --parsable \
    --job-name=liver_6cohort_integration \
    --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/6cohort_integration_%j.out \
    --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/6cohort_integration_%j.err \
    --partition=bigmem \
    --cpus-per-task=32 \
    --mem=400G \
    --time=12:00:00 \
    /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/run_6cohort_integration_slurm.sh)

echo "Integration job submitted: $INTEGRATION_JOB"
echo ""
echo "=== All done: $(date) ==="
