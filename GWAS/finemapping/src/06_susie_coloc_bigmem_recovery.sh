#!/bin/bash -l
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=32
#SBATCH --mem=1500G
#SBATCH --time=90:00:00  # capped at 90h per CLAUDE.md absolute-max policy (T2.3 2026-04-22); resume logic handles re-runs
#SBATCH --job-name=susie_recovery
#SBATCH --output=logs/susie_recovery_%j.out
#SBATCH --error=logs/susie_recovery_%j.err

# Bigmem multi-task SuSiE-COLOC with 3-level resume:
#   1. Task-list level: persistent task file with status per (GWAS, chr)
#   2. Chr level: per-worker gene checkpoints (handled by R script)
#   3. Orchestration: stale worker detection, atomic status, verified merge
#
# Usage:
#   TASK_FILE=data/recovery/tasks_a_safe.txt N_WORKERS=32 sbatch src/06_susie_coloc_bigmem_recovery.sh
#
# Task file format (pipe-separated):
#   GWAS_NAME|CHR|STATUS
#   STATUS in {pending, in_progress, done, failed}

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

RSCRIPT="/gpfs/commons/home/jameslee/micromamba/envs/finemapping/bin/Rscript"

N_WORKERS="${N_WORKERS:-32}"
TASK_FILE="${TASK_FILE:?TASK_FILE not set}"
STATUS_LOCK="${TASK_FILE}.lock"

if [ ! -f "$TASK_FILE" ]; then
  echo "ERROR: TASK_FILE not found: $TASK_FILE"
  exit 1
fi

total_tasks=$(wc -l < "$TASK_FILE")
echo "=== SuSiE-COLOC Recovery: ${TASK_FILE} ==="
echo "Workers per chr: ${N_WORKERS}"
echo "Total tasks: ${total_tasks}"
echo "Job: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo ""

# --- Atomic status update helper ---
update_status() {
  local gwas="$1" chr="$2" new_status="$3"
  flock -x "$STATUS_LOCK" -c "
    awk -F'|' -v g='$gwas' -v c='$chr' -v s='$new_status' \
      'BEGIN{OFS=\"|\"} { if (\$1==g && \$2==c) \$3=s; print }' '$TASK_FILE' > '$TASK_FILE.tmp' && \
    mv '$TASK_FILE.tmp' '$TASK_FILE'
  "
}

# --- Verify SuSiE column presence ---
verify_susie_done() {
  local final="$1"
  [ -f "$final" ] && head -1 "$final" | grep -q "PP.H4.susie"
}

# --- Worker file cleanup ---
cleanup_worker_files() {
  local out_dir="$1" chr="$2"
  rm -f "${out_dir}/worker_chr${chr}_w"*.csv 2>/dev/null || true
  rm -f "${out_dir}/checkpoint_chr${chr}_w"*.csv 2>/dev/null || true
  rm -f "${out_dir}/checkpoint_chr${chr}.csv" 2>/dev/null || true
}

# --- Detect stale worker files (N_WORKERS mismatch) ---
detect_stale_workers() {
  local out_dir="$1" chr="$2" expected_n="$3"
  local existing=$(ls "${out_dir}/worker_chr${chr}_w"*.csv 2>/dev/null | wc -l)
  if [ "$existing" -gt 0 ] && [ "$existing" -ne "$expected_n" ]; then
    echo "  Stale worker files detected (N=$existing, expected=$expected_n) — cleaning"
    cleanup_worker_files "$out_dir" "$chr"
    return 0
  fi
  return 1
}

# --- Main processing loop ---
n_done=0
n_skipped=0
n_failed=0
n_processed=0

while IFS='|' read -r GWAS CHR STATUS; do
  [ -z "$GWAS" ] && continue

  OUT_DIR="results/susie_coloc/${GWAS}"
  OUT_FILE="${OUT_DIR}/susie_coloc_chr${CHR}.csv"
  mkdir -p "$OUT_DIR"

  # Skip if already done
  if [ "$STATUS" = "done" ]; then
    if verify_susie_done "$OUT_FILE"; then
      n_skipped=$((n_skipped + 1))
      continue
    else
      echo "  WARNING: status=done but PP.H4.susie missing — resetting to pending"
      update_status "$GWAS" "$CHR" "pending"
      STATUS="pending"
    fi
  fi

  # Verify column presence even for pending (pre-existing SuSiE-done from other jobs)
  if verify_susie_done "$OUT_FILE"; then
    echo "=== ${GWAS} chr${CHR}: Already SuSiE-complete — marking done ==="
    update_status "$GWAS" "$CHR" "done"
    n_skipped=$((n_skipped + 1))
    continue
  fi

  echo ""
  echo "=== ${GWAS} chr${CHR}: Starting (${N_WORKERS} workers) ==="
  echo "  $(date)"
  echo "  Status on entry: ${STATUS}"

  # Handle stale worker files
  if [ "$STATUS" = "in_progress" ]; then
    echo "  Resuming from prior run"
    detect_stale_workers "$OUT_DIR" "$CHR" "$N_WORKERS" || true
  else
    detect_stale_workers "$OUT_DIR" "$CHR" "$N_WORKERS" || true
  fi

  # Mark in_progress
  update_status "$GWAS" "$CHR" "in_progress"

  # Launch N_WORKERS parallel R processes (stride assignment)
  PIDS=()
  for w in $(seq 0 $((N_WORKERS - 1))); do
    WORKER_ID=${w} N_WORKERS=${N_WORKERS} CHR_FILTER=${CHR} \
      ${RSCRIPT} src/06_susie_coloc.R "${GWAS}" "${CHR}" > "logs/recovery_${SLURM_JOB_ID}_${GWAS}_chr${CHR}_w${w}.out" 2>&1 &
    PIDS+=($!)
  done

  # Wait for all workers
  FAIL=0
  for pid in "${PIDS[@]}"; do
    wait "$pid" || FAIL=$((FAIL + 1))
  done

  if [ $FAIL -eq $N_WORKERS ]; then
    echo "  ERROR: ALL ${N_WORKERS} workers failed for chr${CHR}"
    update_status "$GWAS" "$CHR" "failed"
    n_failed=$((n_failed + 1))
    continue
  elif [ $FAIL -gt 0 ]; then
    echo "  WARNING: ${FAIL}/${N_WORKERS} workers failed (continuing with partial)"
  fi

  # Merge worker outputs
  WORKER_FILES=$(ls ${OUT_DIR}/worker_chr${CHR}_w*.csv 2>/dev/null || true)
  if [ -n "$WORKER_FILES" ]; then
    echo "  Merging $(echo "$WORKER_FILES" | wc -w) worker files..."
    ${RSCRIPT} -e "
      suppressPackageStartupMessages(library(data.table))
      files <- commandArgs(trailingOnly = TRUE)
      dt_list <- lapply(files, function(f) tryCatch(fread(f), error = function(e) NULL))
      dt_list <- dt_list[!sapply(dt_list, is.null)]
      if (length(dt_list) > 0) {
        merged <- rbindlist(dt_list, fill = TRUE)
        merged <- merged[!duplicated(ensembl)]
        fwrite(merged, '${OUT_FILE}')
        cat('  Merged', nrow(merged), 'genes into ${OUT_FILE}\n')
      }
    " ${WORKER_FILES}
  else
    echo "  WARNING: No worker output files found for chr${CHR}"
  fi

  # Verify merge produced SuSiE columns
  if verify_susie_done "$OUT_FILE"; then
    update_status "$GWAS" "$CHR" "done"
    cleanup_worker_files "$OUT_DIR" "$CHR"
    n_done=$((n_done + 1))
    echo "  ${GWAS} chr${CHR}: DONE at $(date)"
  else
    update_status "$GWAS" "$CHR" "failed"
    echo "  ${GWAS} chr${CHR}: merge failed — worker files preserved for retry"
    n_failed=$((n_failed + 1))
  fi

  n_processed=$((n_processed + 1))

done < "$TASK_FILE"

echo ""
echo "=== Summary ==="
echo "  Processed: ${n_processed}"
echo "  Done:      ${n_done}"
echo "  Skipped:   ${n_skipped} (already SuSiE-complete)"
echo "  Failed:    ${n_failed}"
echo "End: $(date)"
