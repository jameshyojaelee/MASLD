#!/bin/bash
# Self-healing watcher for COLOC rerun 19648525 (1,100 tasks).  v2, 2026-08-10.
#
# v1 defect this replaces: it repaired TIMEOUT by resubmitting at the SAME 90h
# wall. That cannot succeed -- a task that needs >90h needs >90h on the retry
# too -- so it would burn another 90h of compute per task while logging
# "repaired". Memory and wall-time are different failures and are now handled
# differently:
#
#   OUT_OF_MEMORY / NODE_FAIL -> resubmit once at 200G (a real repair)
#   TIMEOUT                   -> ALERT only. The wall is already at the 90h
#                                project ceiling, so a resubmit is guaranteed
#                                waste. Needs a human decision.
#
# Also added: a pre-emptive wall warning. chr1 is the risk -- 0 of its 50 tasks
# had completed when this was written, and the eGene-count extrapolation puts it
# near 57h against the 90h wall. Any running task within 8h of the wall is
# flagged with the exact set_timelimit command, which only the USER can run
# (needs sudo), so it surfaces BEFORE the work is lost rather than after.
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
JOB=19648525
SB="$FM/src/seqfunc_v2/finemap/06c_coloc_rerun_full.sbatch"
WALL_H=90
WARN_H=8          # warn when a running task is within this many hours of the wall
declare -A RETRIED
declare -A WARNED
declare -A TIMEDOUT

log() { echo "$(date -u +%FT%TZ) $*"; }

while :; do
  inflight=$(squeue -u "$USER" -h -n susie_coloc -t RUNNING,PENDING 2>/dev/null | wc -l)
  csv=$(find "$FM/results/susie_coloc_rerun" -name 'susie_coloc_chr*.csv' 2>/dev/null | wc -l)

  states=$(sacct -j "$JOB" --format=JobID,State -n -P 2>/dev/null | grep -vE '\.batch|\.extern')

  # --- real repair: memory / node failures -----------------------------------
  oom=$(echo "$states" | awk -F'|' '$2 ~ /OUT_OF_MEMORY|NODE_FAIL/ {split($1,a,"_"); print a[2]}' | sort -u)
  for t in $oom; do
    [ -n "${RETRIED[$t]:-}" ] && continue
    RETRIED[$t]=1
    log "REPAIR task $t (OOM/node-fail) -> resubmit at 200G"
    (cd "$FM" && sbatch --array="$t" --mem=200G --time=${WALL_H}:00:00 "$SB" >/dev/null 2>&1)
  done

  # --- NOT a repair: wall-clock exhaustion -----------------------------------
  tmo=$(echo "$states" | awk -F'|' '$2 ~ /TIMEOUT/ {split($1,a,"_"); print a[2]}' | sort -u)
  for t in $tmo; do
    [ -n "${TIMEDOUT[$t]:-}" ] && continue
    TIMEDOUT[$t]=1
    chr=$(( t % 22 + 1 ))
    log "ALERT task $t (chr${chr}) hit the ${WALL_H}h wall. NOT resubmitting -- the retry"
    log "ALERT   would use the same wall and fail identically. Needs a human decision:"
    log "ALERT   split the chromosome, or raise the wall past the 90h project ceiling."
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
      log "WARN   To rescue it without losing progress, USER runs:"
      log "WARN   sudo /usr/bin/set_timelimit.sh ${JOB}_${tid} 5-00:00:00"
    fi
  done < <(squeue -j "$JOB" -h -t RUNNING -O ArrayTaskID:10,TimeUsed:16 2>/dev/null)

  log "inflight=$inflight csv=$csv/1100 repaired=${#RETRIED[@]} timedout=${#TIMEDOUT[@]}"
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
echo "tasks repaired at 200G: ${#RETRIED[@]}"
echo "tasks lost to the wall: ${#TIMEDOUT[@]}"
echo ""
echo "=== per-chromosome output count (expect 50 each) ==="
for c in $(seq 1 22); do
  n=$(find "$FM/results/susie_coloc_rerun" -name "susie_coloc_chr${c}.csv" 2>/dev/null | wc -l)
  st="OK"; [ "$n" -lt 50 ] && st="*** SHORT"
  printf "  chr%-3s %3d/50 %s\n" "$c" "$n" "$st"
done
