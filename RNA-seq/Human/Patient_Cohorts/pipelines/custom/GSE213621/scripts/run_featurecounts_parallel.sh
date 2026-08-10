#!/bin/bash
#SBATCH --job-name=liver_fc_batch
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/slurm_logs/fc_batch_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/slurm_logs/fc_batch_%A_%a.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --array=0-19

echo "ERROR: legacy canonical featureCounts batches are disabled; use bg001_remediation/recount_array.sbatch." >&2
exit 64

# =============================================================================
# Parallel featureCounts — Array Job (20 batches of ~18 BAMs each)
# =============================================================================
# Each array task processes a subset of BAMs independently.
# A downstream merge script combines all batch outputs into gene_counts.txt.
# =============================================================================

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

# --- Config ---
GTF="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
BAM_LIST="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/scripts/bam_list.txt"
OUTDIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE213621/counts/featurecounts"

TOTAL_BAMS=$(wc -l < "$BAM_LIST")
NUM_BATCHES=20
BATCH_ID=${SLURM_ARRAY_TASK_ID}

# --- Calculate which BAMs this batch processes ---
BATCH_SIZE=$(( (TOTAL_BAMS + NUM_BATCHES - 1) / NUM_BATCHES ))
START_LINE=$(( BATCH_ID * BATCH_SIZE + 1 ))
END_LINE=$(( START_LINE + BATCH_SIZE - 1 ))
if [ "$END_LINE" -gt "$TOTAL_BAMS" ]; then END_LINE=$TOTAL_BAMS; fi

# Extract this batch's BAM paths
BATCH_BAMS=$(sed -n "${START_LINE},${END_LINE}p" "$BAM_LIST")
BATCH_COUNT=$(echo "$BATCH_BAMS" | wc -l)

echo "=== featureCounts Batch $BATCH_ID ==="
echo "BAMs: lines $START_LINE-$END_LINE ($BATCH_COUNT files)"
echo "Started: $(date)"
echo ""

# Skip if no BAMs in this batch (edge case for last batch)
if [ -z "$BATCH_BAMS" ]; then
    echo "No BAMs in this batch. Exiting."
    exit 0
fi

# --- Run featureCounts on this batch ---
mkdir -p "$OUTDIR/batches"
BATCH_OUT="$OUTDIR/batches/batch_${BATCH_ID}_counts.txt"

featureCounts \
    -T ${SLURM_CPUS_PER_TASK} \
    -p --countReadPairs -B \
    -s 2 \
    -a "$GTF" \
    -o "$BATCH_OUT" \
    $BATCH_BAMS

echo ""
echo "Batch $BATCH_ID complete: $(date)"
echo "Output: $BATCH_OUT"
echo "Columns: $(head -2 "$BATCH_OUT" | tail -1 | awk -F'\t' '{print NF}')"
echo "Rows: $(wc -l < "$BATCH_OUT")"
