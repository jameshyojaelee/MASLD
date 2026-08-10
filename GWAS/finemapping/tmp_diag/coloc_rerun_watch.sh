#!/bin/bash
# Self-healing watcher for COLOC rerun 19648525 (1,100 tasks).
#
# Written because the failure modes are KNOWN and untested, not hypothetical:
#   - no LARGE chromosome (chr1 2045 genes, chr2 1284, chr6 925, chr19 1400) has
#     completed yet; everything verified is <=1288 genes
#   - longest completed task is 13:30:47 against a 90h wall; chr1 is 4.4x chr22
#   - the eQTL regeneration OOM'd at 48G on exactly chr1/2/6/11/12
#
# So rather than assume, catch and repair: any task that OOMs or times out is
# resubmitted once at 200G / 90h. Completeness is audited against the expected
# 1,100 outputs, because "array finished" has repeatedly not meant "work done".
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
JOB=19648525
SB="$FM/src/seqfunc_v2/finemap/06c_coloc_rerun_full.sbatch"
declare -A RETRIED

log() { echo "$(date -u +%FT%TZ) $*"; }

while :; do
  inflight=$(squeue -u "$USER" -h -n susie_coloc -t RUNNING,PENDING 2>/dev/null | wc -l)
  csv=$(find "$FM/results/susie_coloc_rerun" -name 'susie_coloc_chr*.csv' 2>/dev/null | wc -l)

  # Repair OOM / TIMEOUT tasks as they appear, without waiting for the array to drain.
  bad=$(sacct -j "$JOB" --format=JobID,State -n -P 2>/dev/null \
        | grep -vE '\.batch|\.extern' \
        | awk -F'|' '$2 ~ /OUT_OF_MEMORY|TIMEOUT|NODE_FAIL/ {split($1,a,"_"); print a[2]}' | sort -u)
  for t in $bad; do
    [ -n "${RETRIED[$t]}" ] && continue
    RETRIED[$t]=1
    log "REPAIR task $t (OOM/timeout) -> resubmit at 200G"
    (cd "$FM" && sbatch --array="$t" --mem=200G --time=90:00:00 "$SB" >/dev/null 2>&1)
  done

  log "inflight=$inflight csv=$csv/1100 repaired=${#RETRIED[@]}"
  [ "$inflight" -eq 0 ] && break
  sleep 900
done

echo ""
echo "=== COLOC RERUN FINISHED ==="
sacct -j "$JOB" --format=State -n 2>/dev/null | grep -vE '\.batch|\.extern' \
  | awk '{print $1}' | sort | uniq -c
echo ""
echo "outputs: $(find "$FM/results/susie_coloc_rerun" -name 'susie_coloc_chr*.csv' 2>/dev/null | wc -l) / 1100"
echo "studies: $(ls -d "$FM"/results/susie_coloc_rerun/*/ 2>/dev/null | wc -l) / 50"
echo "tasks repaired at 200G: ${#RETRIED[@]}"
echo ""
echo "=== per-chromosome output count (expect 50 each) ==="
for c in $(seq 1 22); do
  n=$(find "$FM/results/susie_coloc_rerun" -name "susie_coloc_chr${c}.csv" 2>/dev/null | wc -l)
  st="OK"; [ "$n" -lt 50 ] && st="*** SHORT"
  printf "  chr%-3s %3d/50 %s\n" "$c" "$n" "$st"
done
