#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --job-name=count_r8
#SBATCH --array=1-776%150
#SBATCH --output=logs/alt_contig_recount/count_r8/%a-%A.out
# One library per task from the 776-library manifest; counting itself is run_verified_count.sbatch.
set -euo pipefail
task_root="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas}"
candidate="$task_root/figures/candidates/resource-corrections-20261005-v1"
LIBRARY_MANIFEST="${LIBRARY_MANIFEST:-$candidate/bulk_count_inputs/library_manifest_776.tsv}"
export ANNOTATION_ROOT="${ANNOTATION_ROOT:-$candidate/annotation_resolution_r8}"
COUNT_SAMPLE_ROOT="${COUNT_SAMPLE_ROOT:-$candidate/verified_counts_r8_20261007/samples}"
mkdir -p "$COUNT_SAMPLE_ROOT"
line=$(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1{print $1"\t"$3"\t"$4}' "$LIBRARY_MANIFEST")
IFS=$'\t' read -r SAMPLE_ID SOURCE_BAM layout <<<"$line"
export SAMPLE_ID SOURCE_BAM LIBRARY_LAYOUT="${layout^^}" COUNT_OUTPUT_ROOT="$COUNT_SAMPLE_ROOT/$SAMPLE_ID"
if [[ -s "$COUNT_OUTPUT_ROOT/count_completion.json" ]]; then echo "already counted: $SAMPLE_ID"; exit 0; fi
unset SLURM_TMPDIR
export COUNT_SCRATCH_ROOT=$(mktemp -d "/scratch/count-r8-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}-XXXXXX")
# Task scratch holds two BAM copies per library; remove it when the task ends.
trap 'rm -rf -- "$COUNT_SCRATCH_ROOT"' EXIT
bash "$task_root/scripts/analysis/alt_contig_recount/run_verified_count.sbatch"
