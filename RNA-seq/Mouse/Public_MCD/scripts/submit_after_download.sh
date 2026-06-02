#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 /path/to/dataset_root (e.g., RNA-seq/other_MCD_RNAseq/GSE156918)" >&2
  exit 1
fi

ROOT=$(cd "$1" && pwd)
LOG="$ROOT/download.log"
SUBMIT_LOG="$ROOT/snakemake_submit.log"
DESEQ_LOG="$ROOT/deseq2_submit.log"
MICROMAMBA="/gpfs/commons/groups/sanjana_lab/Cas13/RNA-seq/in-house_MCD_RNAseq/micromamba"
SNAKEMAKE_JOBS="${SNAKEMAKE_JOBS:-200}"
LATENCY_WAIT="${LATENCY_WAIT:-120}"

echo "Waiting for downloads to finish for $ROOT..." >> "$SUBMIT_LOG"
while true; do
  if [[ -f "$LOG" ]] && rg -q "Done\\. Merged FASTQs are in" "$LOG"; then
    if ! pgrep -f "download_and_merge_fastq.sh $ROOT" >/dev/null 2>&1; then
      break
    fi
  fi
  sleep 300
done

echo "Downloads complete. Submitting Snakemake at $(date)" >> "$SUBMIT_LOG"
cd "$ROOT"

$MICROMAMBA run -n rnaseq snakemake \
  -s workflow/Snakefile \
  --jobs "$SNAKEMAKE_JOBS" \
  --latency-wait "$LATENCY_WAIT" \
  --rerun-incomplete \
  --cluster-config workflow/cluster_config.yaml \
  --cluster "sbatch -A {cluster.account} -p {cluster.partition} --time {cluster.time} --mem {cluster.mem} {cluster.extra}" \
  --jobscript workflow/jobscript.sh >> "$SUBMIT_LOG" 2>&1

echo "Submitting DESeq2 at $(date)" >> "$DESEQ_LOG"
sbatch \
  -A nslab \
  -p cpu \
  --qos=nslab \
  --time=06:00:00 \
  --mem=32G \
  -J "deseq2_$(basename "$ROOT")" \
  --output "$ROOT/deseq2_slurm_%j.out" \
  --wrap "cd \"$ROOT\" && $MICROMAMBA run -n rnaseq Rscript scripts/run_deseq2_mcd_only.R" >> "$DESEQ_LOG" 2>&1

echo "Submission complete at $(date)" >> "$SUBMIT_LOG"
