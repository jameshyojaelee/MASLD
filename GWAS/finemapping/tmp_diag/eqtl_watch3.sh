#!/bin/bash
# Self-healing completeness watcher for the eQTL SuSiE regeneration.
#
# Audits against each run's OWN eGene total, parsed from "Processing N eGenes on
# chr C" in its log -- NOT the canonical RDS directory count, which differs by
# ~1-2% and is itself sometimes truncated (canonical chr22 held 404 of 461).
#
# It loops rather than reporting once, because "job COMPLETED" repeatedly turned
# out not to mean "chromosome finished": the 48G array OOM'd on 14 of 22
# chromosomes, and GENE_START_PCT workers only cover their offset-to-end range,
# so a dead base worker silently leaves the head of the gene list unprocessed.
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
MAX_ROUNDS=3

egene_total() {
  grep -ohE "Processing [0-9]+ eGenes on chr $1 " "$FM"/logs/eqtl_susie_*_"$1".out 2>/dev/null \
    | grep -oE "[0-9]+" | head -1
}
find_gaps() {
  local gap="" c cur tot
  for c in $(seq 1 22); do
    cur=$(ls "$FM"/results/eqtl_susie_polyfun/chr${c}/*.rds 2>/dev/null | wc -l)
    tot=$(egene_total "$c"); [ -z "$tot" ] && tot=0
    [ "$tot" -gt 0 ] && [ $((100 * cur / tot)) -lt 95 ] && gap="$gap $c"
  done
  echo "$gap"
}

round=0
while :; do
  n=$(squeue -u "$USER" -h -n eqtl_susie -t RUNNING,PENDING 2>/dev/null | wc -l)
  r=$(find "$FM"/results/eqtl_susie_polyfun -name '*.rds' 2>/dev/null | wc -l)
  echo "$(date -u +%H:%M) in-flight=$n rds=$r"
  if [ "$n" -eq 0 ]; then
    gap=$(find_gaps)
    if [ -z "$gap" ]; then
      echo "ALL CHROMOSOMES >=95% -- REGENERATION COMPLETE"
      break
    fi
    round=$((round + 1))
    if [ "$round" -gt "$MAX_ROUNDS" ]; then
      echo "STILL SHORT after $MAX_ROUNDS auto-resubmits:$gap -- STOPPING, needs a look"
      break
    fi
    echo "auto-resubmit round $round for chromosomes:$gap"
    cd "$FM" || exit 1
    for c in $gap; do
      sbatch --array="$c" --mem=200G \
        --export=ALL,LD_PANEL=polyfun,EQTL_SUSIE_DIR="$FM"/results/eqtl_susie_polyfun \
        src/08_run_eqtl_susie.sh >/dev/null 2>&1
    done
    sleep 90
  fi
  sleep 600
done

echo ""
echo "=== FINAL per-chromosome audit ==="
printf "%-7s %7s %7s %6s %s\n" "chr" "rds" "eGenes" "pct" "STATUS"
short=0
for c in $(seq 1 22); do
  cur=$(ls "$FM"/results/eqtl_susie_polyfun/chr${c}/*.rds 2>/dev/null | wc -l)
  tot=$(egene_total "$c"); [ -z "$tot" ] && tot=0
  pct=0; [ "$tot" -gt 0 ] && pct=$((100 * cur / tot))
  st="OK"
  if [ "$pct" -lt 95 ]; then st="*** SHORT"; short=$((short + 1)); fi
  printf "chr%-4s %7d %7d %5d%% %s\n" "$c" "$cur" "$tot" "$pct" "$st"
done
echo ""
echo "chromosomes below 95%: $short"
echo "TOTAL rds: $(find "$FM"/results/eqtl_susie_polyfun -name '*.rds' 2>/dev/null | wc -l)"
