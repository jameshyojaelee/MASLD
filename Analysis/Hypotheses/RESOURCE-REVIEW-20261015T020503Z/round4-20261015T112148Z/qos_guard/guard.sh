#!/bin/bash
# Moves this round's pending masld_atlas jobs from innovation to nslab; exits (to notify) if one is already running on innovation.
# v2 (2026-10-15 ~15:10 EDT): only jobs whose work directory is under masld_atlas; v1 also matched ~/vc jobs.
LOG=$(dirname $0)/guard.log
while true; do
  squeue -u $USER -h -t PENDING,RUNNING -o "%i|%j|%q|%T|%V|%Z" | while IFS='|' read id name qos state sub wd; do
    [[ "$qos" == innovation && "$name" =~ ^[0-9]{4}$ && "$sub" > "2026-10-15T07:20" && "$wd" == *"/Cas13/masld_atlas"* ]] || continue
    if [[ "$state" == PENDING ]]; then
      if scontrol update JobId=$id Account=nslab QOS=nslab >>$LOG 2>&1; then echo "$(date -Is) moved $id $name to nslab" >>$LOG
      else echo "$(date -Is) FAILED to move $id $name" >>$LOG; echo "STOP $id $name move failed" >>$LOG.flag; fi
    else
      echo "STOP $id $name running on innovation" >>$LOG.flag
    fi
  done
  if [[ -s $LOG.flag ]]; then cat $LOG.flag; exit 0; fi
  sleep 120
done
