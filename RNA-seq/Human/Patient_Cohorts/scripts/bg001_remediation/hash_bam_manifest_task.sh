#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

if [[ $# -ne 2 || "$1" != "--task" || ! "$2" =~ ^[0-7]$ ]]; then
  echo "Usage: $0 --task 0..7" >&2
  exit 64
fi
TASK=$2
: "${BG001_RUN_ROOT:?BG001_RUN_ROOT is required}"
RUN_ROOT=$(realpath "$BG001_RUN_ROOT")
[[ -f "$RUN_ROOT/.bg001_candidate_root" ]] || { echo "Missing candidate sentinel" >&2; exit 64; }
SCRIPT_DIR="$RUN_ROOT/source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
MANIFEST="$RUN_ROOT/manifests/bam_manifest.draft.tsv"
SHARD_DIR="$RUN_ROOT/manifests/bam_hash_shards"
SHARD="$SHARD_DIR/task_${TASK}.tsv"
FAILED="$SHARD_DIR/task_${TASK}.FAILED"
[[ -d "$SHARD_DIR" && ! -e "$SHARD" && ! -e "$FAILED" ]] || {
  echo "Refusing existing or invalid BAM-hash shard destination for task $TASK" >&2
  exit 64
}

python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$RUN_ROOT"

case "$TASK" in
  0) DATASET=GSE130970; EXPECTED_TOTAL=78; START=0; EXPECTED_TASK=78 ;;
  1) DATASET=GSE135251; EXPECTED_TOTAL=216; START=0; EXPECTED_TASK=216 ;;
  2) DATASET=GSE174478; EXPECTED_TOTAL=93; START=0; EXPECTED_TASK=93 ;;
  3) DATASET=GSE240729; EXPECTED_TOTAL=66; START=0; EXPECTED_TASK=66 ;;
  4|5|6|7)
    DATASET=GSE213621
    EXPECTED_TOTAL=367
    CHUNK=$((TASK - 4))
    START=$((CHUNK * 92))
    EXPECTED_TASK=92
    ((CHUNK == 3)) && EXPECTED_TASK=91
    ;;
esac

mapfile -t ALL_ROWS < <(
  awk -F '\t' -v ds="$DATASET" 'NR>1 && $1==ds && $6=="included" {print $2 "\t" $3 "\t" $8 "\t" $9}' "$MANIFEST"
)
[[ ${#ALL_ROWS[@]} -eq $EXPECTED_TOTAL ]] || {
  echo "$DATASET draft manifest expected $EXPECTED_TOTAL included rows, found ${#ALL_ROWS[@]}" >&2
  exit 1
}
ROWS=("${ALL_ROWS[@]:START:EXPECTED_TASK}")
[[ ${#ROWS[@]} -eq $EXPECTED_TASK ]] || { echo "Task $TASK slice cardinality failed" >&2; exit 1; }

TMP="$SHARD.tmp.${SLURM_JOB_ID:-$$}"
trap 'status=$?; if ((status != 0)); then if [[ -f "$TMP" ]]; then mv "$TMP" "$FAILED"; else touch "$FAILED"; fi; fi; exit $status' EXIT
printf 'task\tdataset\tsample_id\tbam_path\tbam_size\tbam_mtime\tbam_sha256\n' > "$TMP"
for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r sample bam expected_size expected_mtime <<< "$row"
  [[ -f "$bam" && ! -L "$bam" ]] || { echo "Missing or symlinked BAM: $bam" >&2; exit 1; }
  size_before=$(stat -c '%s' "$bam")
  mtime_before=$(stat -c '%Y' "$bam")
  [[ "$size_before" == "$expected_size" ]] || { echo "BAM size drift before freeze: $bam" >&2; exit 1; }
  [[ "$mtime_before" == "$expected_mtime" ]] || { echo "BAM mtime drift before freeze: $bam" >&2; exit 1; }
  samtools quickcheck -v "$bam"
  digest=$(sha256sum "$bam" | awk '{print $1}')
  size_after=$(stat -c '%s' "$bam")
  mtime_after=$(stat -c '%Y' "$bam")
  [[ "$size_after" == "$expected_size" && "$mtime_after" == "$expected_mtime" ]] || {
    echo "BAM changed while hashing: $bam" >&2
    exit 1
  }
  [[ "$digest" =~ ^[0-9a-f]{64}$ ]] || { echo "Invalid SHA-256 for $bam" >&2; exit 1; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$TASK" "$DATASET" "$sample" "$bam" "$size_after" "$mtime_after" "$digest" >> "$TMP"
done

python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$RUN_ROOT"
mv "$TMP" "$SHARD"
trap - EXIT
echo "PASS: task $TASK froze ${#ROWS[@]} BAM hashes for $DATASET"
