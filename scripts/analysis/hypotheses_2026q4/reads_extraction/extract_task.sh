#!/bin/bash
# Body shared by run_extract_one.sbatch and run_extract_array.sbatch.
# Usage: extract_task.sh RUN_DIR MANIFEST_INDEX
# The body is one brace group so bash parses the whole file before running it; editing this file
# while tasks are running cannot break a running task.
{
set -euo pipefail
RUN_DIR="${1:?run dir}"
IDX="${2:?manifest index}"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
MANIFEST="$RUN_DIR/library_manifest.tsv"
RES="${RESOURCES:-$RUN_DIR/resources/final}"
# columns: array_index cohort library layout in_pooled_five in_776_counted strandedness_776 bam_bytes bam_path sj_path rpg_path log_final_path
IFS=$'\t' read -r _i COHORT LIB LAYOUT _p _c _s _b BAM SJ RPG LOGF < <(awk -F'\t' -v i="$IDX" 'NR>1 && $1==i {print; exit}' "$MANIFEST")
test -n "${LIB:-}" || { echo "no manifest row for index $IDX" >&2; exit 3; }
LIB_ROOT="${LIB_ROOT:-$RUN_DIR/libraries}"
OUT="$LIB_ROOT/$COHORT/$LIB"
mkdir -p "$LIB_ROOT/$COHORT"
export MASLD_FASTA=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
"$PY" -I "$ROOT/scripts/analysis/hypotheses_2026q4/reads_extraction/extract_reads_one_library.py" \
  --library "$LIB" --cohort "$COHORT" --layout "$LAYOUT" --bam "$BAM" --sj "$SJ" --rpg "$RPG" \
  --log-final "$LOGF" --resources "$RES" --out-dir "$OUT" --threads 2
exit 0
}
