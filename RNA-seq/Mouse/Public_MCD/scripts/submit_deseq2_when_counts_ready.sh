#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 /path/to/dataset_root (e.g., RNA-seq/other_MCD_RNAseq/GSE156918)" >&2
  exit 1
fi

ROOT=$(cd "$1" && pwd)
COUNT="$ROOT/counts/featurecounts/gene_counts.txt"
OUT_TSV="$ROOT/analysis_mcd_vs_control/deseq2_mcd_vs_control.tsv"
LOG="$ROOT/deseq2_submit.log"
MICROMAMBA="/gpfs/commons/groups/sanjana_lab/Cas13/RNA-seq/in-house_MCD_RNAseq/micromamba"
INTERVAL="${INTERVAL:-300}"

echo "Waiting for featureCounts output at $COUNT" >> "$LOG"
while [[ ! -s "$COUNT" ]]; do
  sleep "$INTERVAL"
done

if [[ -s "$OUT_TSV" ]]; then
  echo "DESeq2 output already exists at $OUT_TSV; skipping submission." >> "$LOG"
  exit 0
fi

echo "Counts ready at $(date). Submitting DESeq2." >> "$LOG"
sbatch \
  -A nslab \
  -p cpu \
  --qos=nslab \
  --time=06:00:00 \
  --mem=32G \
  -J "deseq2_$(basename "$ROOT")" \
  --output "$ROOT/deseq2_slurm_%j.out" \
  --wrap "cd \"$ROOT\" && $MICROMAMBA run -n rnaseq Rscript scripts/run_deseq2_mcd_only.R" >> "$LOG" 2>&1
