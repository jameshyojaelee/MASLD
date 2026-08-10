#!/bin/bash
# Self-healing watcher for COLOC rerun 19648525 (1,100 tasks).  v3, 2026-08-10.
#
# CORRECTION TO v2. v2 refused to resubmit a TIMEOUT task, on the reasoning that
# the retry would use the same 90h wall and fail identically. That reasoning was
# wrong because I had not read the script: 06_susie_coloc.R:755 writes
# checkpoint_chr<N>.csv every 200 genes, and :261-273 resumes from it, skipping
# genes already done. 54 such checkpoints are live in the rerun right now. So a
# resubmitted task does NOT start over -- it loses at most the ~5h since its last
# checkpoint and continues. Resubmitting is correct, and v1 was accidentally right.
#
# Failure-specific repair:
#   OUT_OF_MEMORY / NODE_FAIL -> resubmit at 200G (memory is the problem)
#   TIMEOUT                   -> resubmit at the SAME 96G (memory is NOT the
#                                problem; the checkpoint carries the progress).
#                                Up to MAX_TIMEOUT_RETRY rounds, because a task
#                                projected past 90h may need more than one.
#
# Pre-emptive wall warning retained: a rescue via set_timelimit on a RUNNING job
# is still better than a requeue, since it loses nothing at all and avoids
# re-queueing behind the throttle. Only the USER can run it (needs sudo).
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
JOB=19648525
SB="$FM/src/seqfunc_v2/finemap/06c_coloc_rerun_full.sbatch"
WALL_H=90
WARN_H=8
MAX_TIMEOUT_RETRY=3
declare -A RETRIED
declare -A WARNED
declare -A TMO_N

log() { echo "$(date -u +%FT%TZ) $*"; }

while :; do
  inflight=$(squeue -u "$USER" -h -n susie_coloc -t RUNNING,PENDING 2>/dev/null | wc -l)
  csv=$(find "$FM/results/susie_coloc_rerun" -name 'susie_coloc_chr*.csv' 2>/dev/null | wc -l)

  states=$(sacct -j "$JOB" --format=JobID,State -n -P 2>/dev/null | grep -vE '\.batch|\.extern')

  # --- memory / node failure -> more memory ----------------------------------
  for t in $(echo "$states" | awk -F'|' '$2 ~ /OUT_OF_MEMORY|NODE_FAIL/ {split($1,a,"_"); print a[2]}' | sort -u); do
    [ -n "${RETRIED[$t]:-}" ] && continue
    RETRIED[$t]=1
    log "REPAIR task $t (OOM/node-fail) -> resubmit at 200G"
    (cd "$FM" && sbatch --array="$t" --mem=200G --time=${WALL_H}:00:00 "$SB" >/dev/null 2>&1)
  done

  # --- wall exhaustion -> resubmit, checkpoint resumes ------------------------
  for t in $(echo "$states" | awk -F'|' '$2 ~ /TIMEOUT/ {split($1,a,"_"); print a[2]}' | sort -u); do
    n=${TMO_N[$t]:-0}
    [ "$n" -ge "$MAX_TIMEOUT_RETRY" ] && continue
    TMO_N[$t]=$((n+1))
    chr=$(( t % 22 + 1 ))
    log "REPAIR task $t (chr${chr}) hit the ${WALL_H}h wall -> resubmit (round $((n+1))/${MAX_TIMEOUT_RETRY});"
    log "REPAIR   resumes from checkpoint_chr${chr}.csv, losing only since the last 200-gene save"
    (cd "$FM" && sbatch --array="$t" --time=${WALL_H}:00:00 "$SB" >/dev/null 2>&1)
  done

  # --- pre-emptive: running tasks approaching the wall ------------------------
  while read -r tid used; do
    [ -z "$tid" ] && continue
    h=$(echo "$used" | awk -F'[-:]' '{ if (NF==4) print $1*24+$2+$3/60; else if (NF==3) print $1+$2/60; else print $1/60 }')
    over=$(awk -v h="$h" -v w="$WALL_H" -v m="$WARN_H" 'BEGIN{print (h > w-m) ? 1 : 0}')
    if [ "$over" = "1" ] && [ -z "${WARNED[$tid]:-}" ]; then
      WARNED[$tid]=1
      chr=$(( tid % 22 + 1 ))
      log "WARN task $tid (chr${chr}) at ${h%.*}h of ${WALL_H}h wall -- under ${WARN_H}h left."
      log "WARN   Cheapest rescue (loses nothing, no requeue) -- USER runs:"
      log "WARN   sudo /usr/bin/set_timelimit.sh ${JOB}_${tid} 5-00:00:00"
      log "WARN   If not run, this watcher will resubmit on TIMEOUT and resume from checkpoint."
    fi
  done < <(squeue -j "$JOB" -h -t RUNNING -O ArrayTaskID:10,TimeUsed:16 2>/dev/null)

  nt=0; for k in "${!TMO_N[@]}"; do nt=$((nt+1)); done
  log "inflight=$inflight csv=$csv/1100 oom_repaired=${#RETRIED[@]} wall_resubmitted=$nt"
  [ "$inflight" -eq 0 ] && break
  sleep 900
done

echo ""
echo "=== COLOC RERUN FINISHED ==="
sacct -j "$JOB" --format=State -n -P 2>/dev/null | grep -vE '\.batch|\.extern' \
  | awk -F'|' '{print $2}' | sort | uniq -c
echo ""
echo "outputs: $(find "$FM/results/susie_coloc_rerun" -name 'susie_coloc_chr*.csv' 2>/dev/null | wc -l) / 1100"
echo "studies: $(ls -d "$FM"/results/susie_coloc_rerun/*/ 2>/dev/null | wc -l) / 50"
echo "OOM repairs: ${#RETRIED[@]}   wall resubmissions: ${#TMO_N[@]}"
echo ""
echo "=== per-chromosome output count (expect 50 each) ==="
for c in $(seq 1 22); do
  n=$(find "$FM/results/susie_coloc_rerun" -name "susie_coloc_chr${c}.csv" 2>/dev/null | wc -l)
  st="OK"; [ "$n" -lt 50 ] && st="*** SHORT"
  printf "  chr%-3s %3d/50 %s\n" "$c" "$n" "$st"
done
