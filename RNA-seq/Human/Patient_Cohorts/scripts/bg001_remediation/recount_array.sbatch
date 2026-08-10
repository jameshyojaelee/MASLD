#!/usr/bin/env bash
#SBATCH --job-name=featureCounts
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --array=0-7
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
: "${BG001_RUN_ROOT:?BG001_RUN_ROOT is required}"
[[ -f "$BG001_RUN_ROOT/.bg001_candidate_root" ]] || { echo "Missing candidate sentinel" >&2; exit 2; }

SCRIPT_DIR="$BG001_RUN_ROOT/source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
MANIFEST="$BG001_RUN_ROOT/manifests/bam_manifest.tsv"
TASK=${SLURM_ARRAY_TASK_ID:?}

python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$BG001_RUN_ROOT" --require-bam-frozen
[[ -f "$BG001_RUN_ROOT/BASELINE_FROZEN" ]] || { echo "Baseline is not frozen" >&2; exit 64; }

case "$TASK" in
  0) DATASET=GSE130970; OUT="$BG001_RUN_ROOT/counts/$DATASET"; EXTRA=() ;;
  1) DATASET=GSE135251; OUT="$BG001_RUN_ROOT/counts/$DATASET"; EXTRA=() ;;
  2) DATASET=GSE174478; OUT="$BG001_RUN_ROOT/counts/$DATASET"; EXTRA=() ;;
  3) DATASET=GSE240729; OUT="$BG001_RUN_ROOT/counts/$DATASET"; EXTRA=() ;;
  4|5|6|7)
    DATASET=GSE213621
    CHUNK=$((TASK - 4))
    OUT="$BG001_RUN_ROOT/counts/$DATASET/chunks/chunk$CHUNK"
    EXTRA=(--chunk-index "$CHUNK" --chunk-count 4)
    ;;
  *) echo "Unexpected array task: $TASK" >&2; exit 2 ;;
esac

mkdir -p "$BG001_RUN_ROOT/logs"
exec > >(tee "$BG001_RUN_ROOT/logs/recount_${TASK}.out") 2> >(tee "$BG001_RUN_ROOT/logs/recount_${TASK}.err" >&2)
exec bash "$SCRIPT_DIR/run_featurecounts_contract.sh" \
  --dataset "$DATASET" \
  --manifest "$MANIFEST" \
  --output-dir "$OUT" \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  "${EXTRA[@]}"
