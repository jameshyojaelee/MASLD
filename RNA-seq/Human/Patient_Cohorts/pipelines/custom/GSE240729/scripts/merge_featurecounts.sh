#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/logs/fc_merge_%j.log

echo "ERROR: live canonical merge is disabled; use bg001_remediation/merge_validate.sbatch." >&2
exit 64

set -euo pipefail

eval "$(micromamba shell hook -s bash)"
set +u
micromamba activate rnaseq
set -u

OUTDIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE240729/counts/featurecounts"
PERDIR="${OUTDIR}/per_sample"

echo "Merging per-sample featureCounts..."
MERGER="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/shared/scripts/merge_featurecounts_exact.py"
BAM_LIST="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/scripts/bam_list.txt"
SAMPLE_ARGS=()
while IFS= read -r bam; do
    [[ -n "$bam" ]] || continue
    SAMPLE_ARGS+=(--sample-id "$(basename "$(dirname "$bam")")")
done < "$BAM_LIST"
[[ $((${#SAMPLE_ARGS[@]} / 2)) -eq 66 ]] || { echo "Expected 66 explicit samples" >&2; exit 1; }
python3 "$MERGER" \
    --per-sample-dir "$PERDIR" \
    --output "$OUTDIR/gene_counts.txt" \
    --summary "$OUTDIR/gene_counts.txt.summary" \
    --log "$OUTDIR/merge_featurecounts.log" \
    "${SAMPLE_ARGS[@]}"

echo "Done: $(date)"
