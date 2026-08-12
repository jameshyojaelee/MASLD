#!/bin/bash
# Raise the COLOC array throttle to 60 once the verification + profiler jobs are
# done, per the instruction "once the verification is done and we have room for
# more jobs, then increase the throttle to up to 60 jobs".
#
# The throttle was dropped 50 -> 30 to free scheduling room for those three jobs.
# This restores and raises it, automatically, so the overnight window is not lost
# waiting on a human (or on a Claude session) to be awake at the right moment.
#
# "UP TO 60" is read as a ceiling, not a target. 60 tasks x 96 GiB = 5,760 GiB of
# the 7,168 GiB (7 TiB) nslab per-user cap = 80%, leaving ~1,400 GiB for other
# work. If measured headroom at the time cannot support 60, this sets the largest
# safe value instead and says so rather than queueing jobs that would just sit
# PENDING on QOSMaxMemoryPerUser.
JOB=19648525
WATCH_JOBS="19669586 19671233 19720743"
TARGET=60
MEM_PER_TASK_GIB=96
CAP_GIB=7168          # nslab MaxTRESPU mem=7T, confirmed via scontrol show assoc_mgr
RESERVE_GIB=1200      # left free for the user's own jobs
LOG=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/tmp_diag/throttle_raise.log

log() { echo "$(date -u +%FT%TZ) $*" >> "$LOG"; }

log "waiting on verification jobs: $WATCH_JOBS (target throttle $TARGET)"
while :; do
  n=0
  for j in $WATCH_JOBS; do
    n=$(( n + $(squeue -j "$j" -h -t RUNNING,PENDING 2>/dev/null | wc -l) ))
  done
  [ "$n" -eq 0 ] && break
  sleep 300
done
log "verification jobs finished"

# Memory actually committed right now, from SLURM's own accounting rather than a
# guess: MaxTRESPU reports mem=<limit>(<used>) in MB.
used_mb=$(scontrol show assoc_mgr qos=nslab flags=qos 2>/dev/null \
          | grep -o 'MaxTRESPU=cpu=850([0-9]*),mem=[0-9]*([0-9]*)' \
          | head -1 | sed 's/.*mem=[0-9]*(\([0-9]*\)).*/\1/')
[ -z "$used_mb" ] && used_mb=0
used_gib=$(( used_mb / 1024 ))
running=$(squeue -u "$USER" -h -r -n susie_coloc -t RUNNING 2>/dev/null | wc -l)

# How many additional tasks fit inside the cap minus the reserve.
avail_gib=$(( CAP_GIB - RESERVE_GIB - used_gib ))
[ "$avail_gib" -lt 0 ] && avail_gib=0
room=$(( avail_gib / MEM_PER_TASK_GIB ))
safe=$(( running + room ))

if [ "$safe" -ge "$TARGET" ]; then
  new=$TARGET
  log "headroom OK: used=${used_gib}GiB running=${running} room=+${room} -> setting throttle ${new}"
else
  new=$safe
  [ "$new" -lt 30 ] && new=30
  log "headroom LIMITED: used=${used_gib}GiB running=${running} room=+${room};"
  log "  ${TARGET} would exceed the cap minus ${RESERVE_GIB}GiB reserve -> setting throttle ${new} instead"
fi

before=$(scontrol show job "$JOB" 2>/dev/null | grep -o 'ArrayTaskThrottle=[0-9]*' | head -1)
scontrol update JobId="$JOB" ArrayTaskThrottle="$new" 2>>"$LOG"
sleep 5
after=$(scontrol show job "$JOB" 2>/dev/null | grep -o 'ArrayTaskThrottle=[0-9]*' | head -1)
still=$(squeue -u "$USER" -h -r -n susie_coloc -t RUNNING 2>/dev/null | wc -l)
log "throttle ${before} -> ${after}; running ${running} -> ${still} (running tasks are never disturbed by a throttle change)"
log "projected at full throttle: $(( new * MEM_PER_TASK_GIB )) GiB of ${CAP_GIB} GiB cap"
