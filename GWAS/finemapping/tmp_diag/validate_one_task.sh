#!/bin/bash
# Validate the targeted-rerun mechanism on ONE task before committing 1,050.
#
# Three things are asserted, none of which I have actually tested yet:
#   A. seeding a checkpoint makes the script recompute ONLY the excluded genes
#      (if wrong, every task re-does ~900 genes and the "targeted" rerun costs
#       as much as a full one)
#   B. the emitted file still holds every gene, not just the recomputed ones
#      (the merge is rbindlist(old,new)[!duplicated] -- if that misfires the
#       output is truncated and would silently corrupt the results)
#   C. a recovered gene actually flips abf_fallback -> susie
#      (if not, the whole merge bought nothing and the 1,050 tasks are pointless)
#
# It also MEASURES the per-task cost, which so far I have only estimated.
set -euo pipefail
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
cd "$FM"
LOG=tmp_diag/validate_one.log
log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$LOG"; }

while ! grep -q "MERGE + SEED COMPLETE" tmp_diag/merge_rerun.log 2>/dev/null; do
  log "waiting for merge+seed"; sleep 300
done
log "merge complete; validating one task"

# Pick a small chromosome carrying recovered genes: chr21 has exactly 1.
STUDY=$(ls -d results/susie_coloc_rerun/*/ | head -1 | xargs basename)
CHR=21
OUT="results/susie_coloc_rerun/${STUDY}/susie_coloc_chr${CHR}.csv"
CK="results/susie_coloc_rerun/${STUDY}/checkpoint_chr${CHR}.csv"

BEFORE_ROWS=$(( $(wc -l < "$OUT") - 1 ))
CK_ROWS=$(( $(wc -l < "$CK") - 1 ))
EXPECT_RECOMPUTE=$(( BEFORE_ROWS - CK_ROWS ))
log "study=${STUDY} chr=${CHR}"
log "  output rows before      : ${BEFORE_ROWS}"
log "  checkpoint rows (skip)  : ${CK_ROWS}"
log "  genes expected to re-run: ${EXPECT_RECOMPUTE}"
cp -p "$OUT" tmp_diag/validate_before_chr${CHR}.csv

mapfile -t GW < <(awk -F'\t' 'NR>1 && $1!="" {print $1}' config/gwas_registry.tsv | sort -u)
GI=-1; for i in "${!GW[@]}"; do [ "${GW[$i]}" = "$STUDY" ] && GI=$i; done
TASK=$(( GI * 22 + CHR - 1 ))
log "  array task id           : ${TASK}"

J=$(cd "$FM" && sbatch --parsable --array=${TASK} src/seqfunc_v2/finemap/06c_coloc_rerun_full.sbatch)
log "  submitted validation job ${J}"
while squeue -j "$J" -h -t RUNNING,PENDING 2>/dev/null | grep -q .; do sleep 30; done
ELAPSED=$(sacct -j "$J" --format=ElapsedRaw -n -P 2>/dev/null | head -1)
log "  finished in ${ELAPSED}s"

LG=$(ls -t logs/coloc_rerun/coloc_${J}_${TASK}.out 2>/dev/null | head -1)
RESUMED=$(grep -m1 -o "Resuming from checkpoint: [0-9]* genes" "$LG" 2>/dev/null || echo "NONE")
TESTED=$(grep -o "Genes tested: [0-9]*" "$LG" 2>/dev/null | tail -1)
AFTER_ROWS=$(( $(wc -l < "$OUT") - 1 ))

echo "" | tee -a "$LOG"
echo "=================== VALIDATION ===================" | tee -a "$LOG"
echo "  A. only the excluded genes recomputed?" | tee -a "$LOG"
echo "     $RESUMED  (expected ${CK_ROWS})" | tee -a "$LOG"
echo "  B. output still complete?" | tee -a "$LOG"
echo "     rows ${BEFORE_ROWS} -> ${AFTER_ROWS}   $( [ "$BEFORE_ROWS" -eq "$AFTER_ROWS" ] && echo OK || echo '*** ROW COUNT CHANGED' )" | tee -a "$LOG"
echo "  C. did a recovered gene flip to susie?" | tee -a "$LOG"
python3 - "$OUT" tmp_diag/validate_before_chr${CHR}.csv tmp_diag/recovered_genes.tsv "$CHR" <<'PY' 2>&1 | tee -a "$LOG"
import csv,sys
new,old,rec,chrom=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
want={l.split("\t")[1].strip() for l in open(rec) if l.split("\t")[0]==f"chr{chrom}"}
def load(p): return {r["ensembl"]:r for r in csv.DictReader(open(p))}
a,b=load(old),load(new)
for g in sorted(want):
    if g in a and g in b:
        print(f"     {g}  method {a[g]['method']} -> {b[g]['method']}   "
              f"PP.H4.susie {a[g]['PP.H4.susie'] or 'NA'} -> {b[g]['PP.H4.susie'] or 'NA'}")
PY
echo "  D. measured cost: ${ELAPSED}s/task  ->  1050 tasks = $(( ELAPSED * 1050 / 3600 )) task-hours" | tee -a "$LOG"
echo "==================================================" | tee -a "$LOG"
