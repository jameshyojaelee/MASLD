#!/bin/bash
# Self-healing watcher for COLOC rerun 19648525 (1,100 tasks).  v4, 2026-08-13.
#
# CORRECTION TO v3. v3 repaired only OUT_OF_MEMORY / NODE_FAIL / TIMEOUT. It
# therefore sat silent through 41 FAILED tasks:
#   * 40 from a ~35-min window on 2026-08-11 when 06_susie_coloc.R was
#     transiently corrupted (file.remove -> ve). Those died on the LAST line,
#     AFTER fwrite, so all 40 outputs are complete and verified well-formed.
#     Re-running them would burn ~166 task-hours to regenerate identical files.
#   * 1 (task 915, PanUKBB_AFR_ALT chr14) was a Bus error inside eigen() on an
#     otherwise-healthy node -- a transient hardware fault that killed the run
#     mid-flight, leaving NO output.
#
# Those two need opposite responses, and no SLURM state distinguishes them: both
# are "FAILED", exit 1:0 and 7:0. What distinguishes them is whether the WORK
# LANDED. So repair is now decided on the OUTPUT, not the state:
#
#     any task with no output file  -> resubmit (checkpoint resumes it)
#     any task with an output file  -> leave alone, whatever SLURM calls it
#
# That is strictly better than state-matching: it also catches states nobody
# thought to enumerate, and it can never re-run work that already succeeded.
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
JOB=19648525
SB="$FM/src/seqfunc_v2/finemap/06c_coloc_rerun_full.sbatch"
WALL_H=90
MAX_RETRY=3
declare -A RETRY
declare -A OOM_SEEN

log() { echo "$(date -u +%FT%TZ) $*"; }

# Same study list the sbatch derives, so task -> output path matches exactly.
mapfile -t GWAS_LIST < <(awk -F'\t' 'NR>1 && $1!="" {print $1}' "${FM}/config/gwas_registry.tsv" | sort -u)

out_path_for() {   # $1 = array task id
  local t=$1 gi=$(( $1 / 22 )) chr=$(( $1 % 22 + 1 ))
  [ "$gi" -ge "${#GWAS_LIST[@]}" ] && { echo ""; return; }
  echo "${FM}/results/susie_coloc_rerun/${GWAS_LIST[$gi]}/susie_coloc_chr${chr}.csv"
}

while :; do
  inflight=$(squeue -u "$USER" -h -r -n susie_coloc -t RUNNING,PENDING 2>/dev/null | wc -l)
  csv=$(find "$FM/results/susie_coloc_rerun" -name 'susie_coloc_chr*.csv' 2>/dev/null | wc -l)

  states=$(sacct -j "$JOB" --format=JobID,State -n -P 2>/dev/null | grep -vE '\.batch|\.extern')

  # Any terminal-but-not-COMPLETED task. Deliberately broad -- the output check
  # below is what decides, so an over-broad match here is harmless.
  bad=$(echo "$states" | awk -F'|' '$2 !~ /COMPLETED|RUNNING|PENDING|REQUEUED|RESIZING|SUSPENDED/ {split($1,a,"_"); print a[2]}' | sort -un)

  for t in $bad; do
    [ -z "$t" ] && continue
    # Already queued or running again (e.g. an earlier repair)? leave it.
    squeue -u "$USER" -h -r -t RUNNING,PENDING -O ArrayTaskID:10 2>/dev/null | tr -d ' ' | grep -qx "$t" && continue
    op=$(out_path_for "$t"); [ -z "$op" ] && continue
    if [ -f "$op" ]; then
      continue                       # work landed -- never re-run it
    fi
    n=${RETRY[$t]:-0}
    if [ "$n" -ge "$MAX_RETRY" ]; then
      [ -z "${OOM_SEEN[$t]:-}" ] && { OOM_SEEN[$t]=1; log "GIVE UP task $t after ${MAX_RETRY} retries -- needs a human"; }
      continue
    fi
    RETRY[$t]=$((n+1))
    # Memory failures get more memory; everything else keeps 96G and relies on
    # the 200-gene checkpoint to resume where it died.
    if echo "$states" | grep -q "_${t}|OUT_OF_MEMORY"; then
      log "REPAIR task $t (OOM, no output) -> resubmit at 200G (try $((n+1))/${MAX_RETRY})"
      (cd "$FM" && sbatch --array="$t" --mem=200G --time=${WALL_H}:00:00 "$SB" >/dev/null 2>&1)
    else
      log "REPAIR task $t (terminal, no output) -> resubmit (try $((n+1))/${MAX_RETRY}); resumes from checkpoint"
      (cd "$FM" && sbatch --array="$t" --time=${WALL_H}:00:00 "$SB" >/dev/null 2>&1)
    fi
  done

  nr=0; for k in "${!RETRY[@]}"; do nr=$((nr+1)); done
  log "inflight=$inflight csv=$csv/1100 repaired=$nr"
  [ "$inflight" -eq 0 ] && [ "$csv" -ge 1100 ] && break
  # Do not exit on an empty queue while outputs are still missing -- that is
  # exactly the state a repair is for. Give the repair a cycle to appear.
  if [ "$inflight" -eq 0 ] && [ "$csv" -lt 1100 ]; then
    log "queue empty but only $csv/1100 outputs -- holding for repair"
    sleep 300
    continue
  fi
  sleep 900
done

echo ""
echo "=== COLOC RERUN FINISHED ==="
echo "outputs: $(find "$FM/results/susie_coloc_rerun" -name 'susie_coloc_chr*.csv' 2>/dev/null | wc -l) / 1100"
echo "studies: $(ls -d "$FM"/results/susie_coloc_rerun/*/ 2>/dev/null | wc -l) / 50"
echo "repairs: ${#RETRY[@]}"
for c in $(seq 1 22); do
  n=$(find "$FM/results/susie_coloc_rerun" -name "susie_coloc_chr${c}.csv" 2>/dev/null | wc -l)
  st="OK"; [ "$n" -lt 50 ] && st="*** SHORT"
  printf "  chr%-3s %3d/50 %s\n" "$c" "$n" "$st"
done
