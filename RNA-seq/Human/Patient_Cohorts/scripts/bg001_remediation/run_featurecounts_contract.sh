#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

usage() {
  echo "Usage: $0 --dataset ID --manifest TSV --output-dir DIR [--chunk-index N --chunk-count 4] [--threads 8]" >&2
  exit 2
}

DATASET=""
MANIFEST=""
OUTPUT_DIR=""
CHUNK_INDEX=""
CHUNK_COUNT=""
THREADS=8
EXPECTED_GTF_SHA="73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
GTF="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"

while (($#)); do
  case "$1" in
    --dataset) DATASET=${2:?}; shift 2 ;;
    --manifest) MANIFEST=${2:?}; shift 2 ;;
    --output-dir) OUTPUT_DIR=${2:?}; shift 2 ;;
    --chunk-index) CHUNK_INDEX=${2:?}; shift 2 ;;
    --chunk-count) CHUNK_COUNT=${2:?}; shift 2 ;;
    --threads) THREADS=${2:?}; shift 2 ;;
    *) usage ;;
  esac
done

[[ -n "$DATASET" && -n "$MANIFEST" && -n "$OUTPUT_DIR" ]] || usage
[[ "$OUTPUT_DIR" = /* && "$MANIFEST" = /* ]] || { echo "Manifest and output paths must be absolute" >&2; exit 2; }
[[ -f "$MANIFEST" ]] || { echo "Missing manifest: $MANIFEST" >&2; exit 2; }
[[ ! -e "$OUTPUT_DIR" ]] || { echo "Refusing existing output directory: $OUTPUT_DIR" >&2; exit 2; }
RUN_ROOT=$(realpath "$(dirname "$MANIFEST")/..")
[[ -f "$RUN_ROOT/.bg001_candidate_root" && ! -L "$RUN_ROOT/.bg001_candidate_root" ]] || {
  echo "Manifest is not inside a marked BG-001 candidate root" >&2
  exit 2
}
[[ $(realpath "$MANIFEST") == "$RUN_ROOT/manifests/bam_manifest.tsv" ]] || {
  echo "Execution requires the candidate's exact frozen BAM manifest" >&2
  exit 2
}
[[ -f "$RUN_ROOT/BAM_MANIFEST_FROZEN" ]] || { echo "BAM manifest is not frozen" >&2; exit 2; }
OUTPUT_RESOLVED=$(realpath -m "$OUTPUT_DIR")
case "$OUTPUT_RESOLVED" in
  "$RUN_ROOT"/counts/*) ;;
  *) echo "Candidate recount output escapes RUN_ROOT/counts" >&2; exit 2 ;;
esac
[[ ! -e "$OUTPUT_RESOLVED" ]] || { echo "Refusing existing resolved output directory: $OUTPUT_RESOLVED" >&2; exit 2; }

case "$DATASET" in
  GSE130970) EXPECTED_TOTAL=78 ;;
  GSE135251) EXPECTED_TOTAL=216 ;;
  GSE174478) EXPECTED_TOTAL=93 ;;
  GSE213621) EXPECTED_TOTAL=367 ;;
  GSE240729) EXPECTED_TOTAL=66 ;;
  *) echo "Dataset is outside the BG-001 recount contract: $DATASET" >&2; exit 2 ;;
esac

if [[ -n "$CHUNK_INDEX" || -n "$CHUNK_COUNT" ]]; then
  [[ "$DATASET" == GSE213621 && "$CHUNK_COUNT" == 4 && "$CHUNK_INDEX" =~ ^[0-3]$ ]] || {
    echo "Only GSE213621 supports --chunk-index 0..3 --chunk-count 4" >&2
    exit 2
  }
fi

eval "$(micromamba shell hook --shell=bash)"
set +u
micromamba activate rnaseq
set -u

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$RUN_ROOT" --require-bam-frozen
FC_VERSION=$({ featureCounts -v 2>&1 || true; } | tr -d '\r' | awk 'NF {last=$0} END {print last}')
[[ "$FC_VERSION" == "featureCounts v2.1.1" ]] || { echo "Expected exact featureCounts v2.1.1, observed: $FC_VERSION" >&2; exit 1; }
GTF_SHA=$(sha256sum "$GTF" | awk '{print $1}')
[[ "$GTF_SHA" == "$EXPECTED_GTF_SHA" ]] || { echo "GENCODE v49 GTF hash mismatch" >&2; exit 1; }

mkdir -p "$OUTPUT_RESOLVED"
OUTPUT_DIR="$OUTPUT_RESOLVED"
trap 'status=$?; if ((status != 0)); then touch "$OUTPUT_DIR/FAILED"; fi; exit $status' EXIT

FC_EXE=$(realpath "$(command -v featureCounts)")
SAMTOOLS_EXE=$(realpath "$(command -v samtools)")
FC_EXE_SHA=$(sha256sum "$FC_EXE" | awk '{print $1}')
SAMTOOLS_EXE_SHA=$(sha256sum "$SAMTOOLS_EXE" | awk '{print $1}')
SAMTOOLS_VERSION=$(samtools --version | awk 'NR==1 {print; exit}')
ENVIRONMENT="$OUTPUT_DIR/environment.txt"
{
  printf 'featureCounts_version\t%s\n' "$FC_VERSION"
  printf 'featureCounts_executable\t%s\n' "$FC_EXE"
  printf 'featureCounts_executable_sha256\t%s\n' "$FC_EXE_SHA"
  printf 'samtools_version\t%s\n' "$SAMTOOLS_VERSION"
  printf 'samtools_executable\t%s\n' "$SAMTOOLS_EXE"
  printf 'samtools_executable_sha256\t%s\n' "$SAMTOOLS_EXE_SHA"
  printf 'gtf_path\t%s\n' "$GTF"
  printf 'gtf_sha256\t%s\n' "$GTF_SHA"
  printf 'micromamba_explicit_begin\n'
  micromamba list --explicit --sha256
  printf 'micromamba_explicit_end\n'
} > "$ENVIRONMENT"

mapfile -t ALL_ROWS < <(awk -F '\t' -v ds="$DATASET" 'NR>1 && $1==ds && $6=="included" {print $2 "\t" $3 "\t" $8 "\t" $9 "\t" $10}' "$MANIFEST")
[[ ${#ALL_ROWS[@]} -eq $EXPECTED_TOTAL ]] || {
  echo "$DATASET manifest expected $EXPECTED_TOTAL included rows, found ${#ALL_ROWS[@]}" >&2
  exit 1
}

START=0
END=$EXPECTED_TOTAL
if [[ -n "$CHUNK_INDEX" ]]; then
  START=$((CHUNK_INDEX * 92))
  END=$((START + 92))
  ((END > EXPECTED_TOTAL)) && END=$EXPECTED_TOTAL
fi
ROWS=("${ALL_ROWS[@]:START:END-START}")
EXPECTED_OUTPUT=${#ROWS[@]}

BAMS=()
: > "$OUTPUT_DIR/bam_checksums.tsv"
printf 'dataset\tsample_id\tbam_path\tbam_size\tbam_mtime\tbam_sha256\n' > "$OUTPUT_DIR/bam_checksums.tsv"
for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r sample bam expected_size expected_mtime expected_sha <<< "$row"
  expected_sha=${expected_sha%$'\r'}
  [[ -f "$bam" ]] || { echo "Missing BAM: $bam" >&2; exit 1; }
  size=$(stat -c '%s' "$bam")
  mtime=$(stat -c '%Y' "$bam")
  [[ "$size" == "$expected_size" ]] || { echo "BAM size drift: $bam" >&2; exit 1; }
  [[ "$mtime" == "$expected_mtime" ]] || { echo "BAM mtime drift: $bam" >&2; exit 1; }
  samtools quickcheck -v "$bam"
  sha=$(sha256sum "$bam" | awk '{print $1}')
  [[ "$expected_sha" =~ ^[0-9a-f]{64}$ && "$sha" == "$expected_sha" ]] || {
    echo "BAM SHA-256 drift after manifest freeze: $bam" >&2
    exit 1
  }
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$DATASET" "$sample" "$bam" "$size" "$mtime" "$sha" >> "$OUTPUT_DIR/bam_checksums.tsv"
  BAMS+=("$bam")
done

COUNTS="$OUTPUT_DIR/gene_counts.txt"
LOG="$OUTPUT_DIR/featureCounts.log"
featureCounts -T "$THREADS" -p --countReadPairs -B -s 2 -a "$GTF" -o "$COUNTS" "${BAMS[@]}" >"$LOG" 2>&1

for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r sample bam expected_size expected_mtime expected_sha <<< "$row"
  expected_sha=${expected_sha%$'\r'}
  [[ $(stat -c '%s' "$bam") == "$expected_size" && $(stat -c '%Y' "$bam") == "$expected_mtime" ]] || {
    echo "BAM changed during recount: $bam" >&2
    exit 1
  }
  sha_after=$(sha256sum "$bam" | awk '{print $1}')
  [[ "$sha_after" == "$expected_sha" ]] || {
    echo "BAM SHA-256 changed during recount: $bam" >&2
    exit 1
  }
done

python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$RUN_ROOT" --require-bam-frozen
VALIDATE_EXTRA=()
[[ -n "$CHUNK_INDEX" ]] && VALIDATE_EXTRA=(--chunk-index "$CHUNK_INDEX")
python3 "$SCRIPT_DIR/validate_featurecounts.py" \
  --dataset "$DATASET" \
  --manifest "$MANIFEST" \
  --counts "$COUNTS" \
  --summary "$COUNTS.summary" \
  --log "$LOG" \
  --environment "$ENVIRONMENT" \
  --expected-samples "$EXPECTED_OUTPUT" \
  --output "$OUTPUT_DIR/validation.json" \
  "${VALIDATE_EXTRA[@]}"

sha256sum \
  "$OUTPUT_DIR/bam_checksums.tsv" \
  "$OUTPUT_DIR/environment.txt" \
  "$OUTPUT_DIR/gene_counts.txt" \
  "$OUTPUT_DIR/gene_counts.txt.summary" \
  "$OUTPUT_DIR/featureCounts.log" \
  "$OUTPUT_DIR/validation.json" \
  > "$OUTPUT_DIR/provenance.sha256"
touch "$OUTPUT_DIR/COMPLETE"
trap - EXIT
