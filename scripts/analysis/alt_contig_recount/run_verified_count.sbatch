#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# One explicitly selected library; no array. Original BAM and named counts remain read-only.
set -euo pipefail
: "${SOURCE_BAM:?}" "${SAMPLE_ID:?}" "${LIBRARY_LAYOUT:?}" "${ANNOTATION_ROOT:?}" "${COUNT_OUTPUT_ROOT:?}"
task_root="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas}"
binary_root=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
"$binary_root/python" -B - "$ANNOTATION_ROOT/resolution_summary.json" <<'PY'
import json,sys
s=json.load(open(sys.argv[1]))
assert s.get('annotation_resolution_complete') is True and s.get('n_unresolved_annotations') == 0, 'Full annotation correspondence unresolved'
PY
test -s "$ANNOTATION_ROOT/grouped_counting_annotation.gtf.gz"
test -s "$SOURCE_BAM"
test ! -e "$COUNT_OUTPUT_ROOT"
[[ "$LIBRARY_LAYOUT" == SINGLE || "$LIBRARY_LAYOUT" == PAIRED ]]
scratch_parent="${SLURM_TMPDIR:-${COUNT_SCRATCH_ROOT:?COUNT_SCRATCH_ROOT required without SLURM_TMPDIR}}"
task_scratch=$(mktemp -d "$scratch_parent/resource-count-${SAMPLE_ID}-XXXXXXXX")
flags=(-s 2 -M -R BAM)
if [[ "$LIBRARY_LAYOUT" == PAIRED ]]; then flags+=(-p --countReadPairs -B); fi
"$binary_root/featureCounts" -T "${SLURM_CPUS_PER_TASK:-4}" "${flags[@]}" \
  -a "$ANNOTATION_ROOT/grouped_counting_annotation.gtf.gz" \
  -o "$task_scratch/alignment_counts_not_for_analysis.tsv" --Rpath "$task_scratch" "$SOURCE_BAM"
assigned="$task_scratch/$(basename "$SOURCE_BAM").featureCounts.bam"
test -s "$assigned"
"$binary_root/samtools" sort -N -@ 2 -m 1G -T "$task_scratch/name-sort" \
  -o "$task_scratch/name-sorted.bam" "$assigned"
"$binary_root/python" -B "$task_root/scripts/analysis/alt_contig_recount/count_verified_fragments.py" \
  --annotation-summary "$ANNOTATION_ROOT/resolution_summary.json" \
  --assignment-bam "$task_scratch/name-sorted.bam" --sample-id "$SAMPLE_ID" \
  --layout "$LIBRARY_LAYOUT" --output "$COUNT_OUTPUT_ROOT"
# Scratch is retained for inspection; scheduler-owned scratch follows its normal lifecycle.
