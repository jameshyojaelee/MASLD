FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
# Wait for ALL eqtl_susie work to clear, then audit completeness per chromosome
# against the canonical gene set -- a job COMPLETING is not the same as the
# chromosome being finished, which is how the chr7/11/12/20 OOM gaps were missed.
while :; do
  n=$(squeue -u $USER -h -n eqtl_susie -t RUNNING,PENDING 2>/dev/null | wc -l)
  r=$(find $FM/results/eqtl_susie_polyfun -name '*.rds' 2>/dev/null | wc -l)
  echo "$(date -u +%H:%M) in-flight=$n rds=$r/18877 ($((100*r/18877))%)"
  [ "$n" -eq 0 ] && break
  sleep 600
done
echo ""
echo "=== ALL EQTL JOBS CLEARED — completeness audit ==="
printf "%-7s %7s %7s %6s %s\n" "chr" "new" "canon" "pct" "STATUS"
short=0
for c in $(seq 1 22); do
  cur=$(ls $FM/results/eqtl_susie_polyfun/chr${c}/*.rds 2>/dev/null|wc -l)
  can=$(ls $FM/results/eqtl_susie/chr${c}/*.rds 2>/dev/null|wc -l)
  pct=0; [ "$can" -gt 0 ] && pct=$((100*cur/can))
  st="OK"; if [ "$pct" -lt 95 ]; then st="*** SHORT"; short=$((short+1)); fi
  printf "chr%-4s %7d %7d %5d%% %s\n" "$c" "$cur" "$can" "$pct" "$st"
done
echo ""
echo "chromosomes below 95%: $short"
echo "=== OOM / FAILED across every eqtl job this session ==="
for j in 19559148 19567542 19568514 19569923 19570061 19570062 19570063 19570064 \
         19570065 19570066 19570067 19570068 19570069 19570070 19570071 19570072 \
         19570073 19570074 19570075 19570076 19570077 19595717 19595718 19595719 \
         19595720 19595721; do
  sacct -j $j --format=JobID,State -n -P 2>/dev/null | grep -vE "\.batch|\.extern"
done | awk -F'|' '$2 ~ /OUT_OF|FAILED|TIMEOUT|CANCELLED/ {print "  " $1, $2}' | sort -u
echo "(empty above = no failures)"
