#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 /path/to/dataset_root (e.g., RNA-seq/other_MCD_RNAseq/GSE156918)" >&2
  exit 1
fi

ROOT=$(cd "$1" && pwd)
META="$ROOT/metadata/samples.tsv"
RUNS="$ROOT/metadata/sra_runs.txt"
RAW_DIR="$ROOT/fastq/raw"
MERGED_DIR="$ROOT/fastq"
PARALLEL_DOWNLOADS="${PARALLEL_DOWNLOADS:-4}"

if [[ ! -f "$META" ]]; then
  echo "Missing $META" >&2
  exit 1
fi
if [[ ! -f "$RUNS" ]]; then
  echo "Missing $RUNS" >&2
  exit 1
fi

mkdir -p "$RAW_DIR" "$MERGED_DIR"

echo "Downloading SRA runs listed in $RUNS (parallel=${PARALLEL_DOWNLOADS})..."
tmp_list=$(mktemp)
cleanup() {
  rm -f "$tmp_list"
}
trap cleanup EXIT
while read -r run; do
  [[ -z "$run" ]] && continue
  info=$(curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n +2)
  if [[ -z "$info" ]]; then
    echo "No ENA entry for $run" >&2
    continue
  fi
  ftp_paths=$(echo "$info" | cut -f2)
  IFS=';' read -r -a files <<< "$ftp_paths"
  for f in "${files[@]}"; do
    [[ -z "$f" ]] && continue
    fname=$(basename "$f")
    dest="$RAW_DIR/$fname"
    if [[ ! -s "$dest" ]]; then
      echo "https://$f $dest" >> "$tmp_list"
    fi
  done
done < "$RUNS"

if [[ -s "$tmp_list" ]]; then
  echo "Queueing downloads..."
  active=0
  while read -r url dest; do
    [[ -z "$url" ]] && continue
    fname=$(basename "$url")
    echo "  downloading $fname"
    (
      curl -L --retry 5 --retry-delay 5 --retry-all-errors "$url" -o "$dest"
    ) &
    active=$((active + 1))
    if [[ "$active" -ge "$PARALLEL_DOWNLOADS" ]]; then
      wait -n
      active=$((active - 1))
    fi
  done < "$tmp_list"
  wait
else
  echo "All FASTQs already present in $RAW_DIR"
fi

echo "Merging runs into per-sample FASTQs..."
tail -n +2 "$META" | while IFS=$'\t' read -r sample_id diet diet_detail genotype biosample runs layout fastq_r1 fastq_r2 title; do
  [[ -z "$sample_id" ]] && continue
  IFS=';' read -r -a run_list <<< "$runs"
  if [[ "$layout" == "PAIRED" ]]; then
    out1="$MERGED_DIR/${sample_id}_1.fastq.gz"
    out2="$MERGED_DIR/${sample_id}_2.fastq.gz"
    if [[ ! -s "$out1" ]]; then
      files1=()
      for run in "${run_list[@]}"; do
        files1+=("$RAW_DIR/${run}_1.fastq.gz")
      done
      cat "${files1[@]}" > "$out1"
    fi
    if [[ ! -s "$out2" ]]; then
      files2=()
      for run in "${run_list[@]}"; do
        files2+=("$RAW_DIR/${run}_2.fastq.gz")
      done
      cat "${files2[@]}" > "$out2"
    fi
  else
    out1="$MERGED_DIR/${sample_id}.fastq.gz"
    if [[ ! -s "$out1" ]]; then
      files1=()
      for run in "${run_list[@]}"; do
        files1+=("$RAW_DIR/${run}.fastq.gz")
      done
      cat "${files1[@]}" > "$out1"
    fi
  fi
done

echo "Done. Merged FASTQs are in $MERGED_DIR"
