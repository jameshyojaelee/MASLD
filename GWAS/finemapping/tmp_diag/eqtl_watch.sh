FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
while :; do
  run=$(squeue -u $USER -h -n eqtl_susie -t RUNNING,PENDING 2>/dev/null | wc -l)
  n=$(find $FM/results/eqtl_susie_polyfun -name '*.rds' 2>/dev/null | wc -l)
  echo "$(date -u +%H:%M) eqtl_susie in-flight=$run  rds=$n / 18877 ($((100*n/18877))%)"
  [ "$run" -eq 0 ] && { echo "ALL EQTL TASKS FINISHED"; break; }
  sleep 600
done
echo "=== final per-chr ==="
for c in $(seq 1 22); do printf "chr%-3s %5d  " "$c" "$(ls $FM/results/eqtl_susie_polyfun/chr${c}/*.rds 2>/dev/null|wc -l)"; [ $((c%5)) -eq 0 ] && echo ""; done; echo ""
echo "=== OOM/FAILED across all eqtl jobs ==="
for j in 19559148 19567542 19568514 19569923; do sacct -j $j --format=JobID,State -n 2>/dev/null | grep -vE "\.batch|\.extern" | grep -E "OUT_OF|FAILED|TIMEOUT"; done
