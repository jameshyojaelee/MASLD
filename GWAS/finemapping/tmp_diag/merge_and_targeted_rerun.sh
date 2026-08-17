#!/bin/bash
# Merge the recovered eQTL fits into production, then re-run COLOC for ONLY the
# gene-study pairs that can possibly change.  2026-08-16.
#
# WHY A TARGETED RERUN IS POSSIBLE AT ALL
# 06_susie_coloc.R:345 skips any gene listed in checkpoint_chr<N>.csv, and the
# final merge (:817-826) does rbindlist(old, new)[!duplicated(ensembl)] -- OLD
# rows win. So seeding the checkpoint with the existing output MINUS the
# recovered genes makes the script recompute exactly those genes and nothing
# else, then re-emit the complete file. No new analysis code; the script's own
# resume machinery does it.
#
# COST: 144 genes x 50 studies = 7,200 gene-computations against 948,503 for a
# full rerun (0.76%). Per-task startup (GWAS + eQTL load) dominates, so this is
# ~1,100 short tasks rather than ~6,000 task-hours.
#
# SAFETY
#   * production fits for all 204 genes are BACKED UP before anything is copied
#   * only the 144 CONVERGED staged fits are merged; the 60 failures are left
#     exactly as they were
#   * the rerun writes to results/susie_coloc_rerun/ -- canonical is untouched
#   * a pre-merge inventory is recorded so the merge is reversible
set -euo pipefail
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
cd "$FM"
LOG=tmp_diag/merge_rerun.log
log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$LOG"; }

BACKUP="results/eqtl_susie_polyfun_backup_2026-08-16"
REC=tmp_diag/recovered_genes.tsv        # <chrN>\t<ENSG>, converged fits only

# ---- 0. wait for the refit array to finish -------------------------------
while squeue -u "$USER" -h -r -n susieR -t RUNNING,PENDING 2>/dev/null | grep -q .; do
  log "waiting for the refit array to drain"; sleep 300
done

n_rec=$(wc -l < "$REC")
log "recovered fits to merge: ${n_rec}"

# ---- 1. back up the production fits we are about to overwrite ------------
mkdir -p "$BACKUP"
n_bk=0
while IFS=$'\t' read -r d g; do
  src="results/eqtl_susie_polyfun/${d}/${g}_susie.rds"
  if [ -f "$src" ]; then mkdir -p "$BACKUP/$d"; cp -p "$src" "$BACKUP/$d/"; n_bk=$((n_bk+1)); fi
done < "$REC"
log "backed up ${n_bk} production fits -> ${BACKUP}"

# ---- 2. merge ONLY the converged recovered fits --------------------------
n_mg=0
while IFS=$'\t' read -r d g; do
  src="results/eqtl_susie_refit/${d}/${g}_susie.rds"
  dst="results/eqtl_susie_polyfun/${d}/${g}_susie.rds"
  if [ -f "$src" ]; then cp -p "$src" "$dst"; n_mg=$((n_mg+1)); fi
done < "$REC"
log "merged ${n_mg} recovered fits into production"
log "production fit count now: $(find results/eqtl_susie_polyfun -name '*_susie.rds' | wc -l)"

# ---- 3. seed checkpoints so only the recovered genes recompute -----------
# For every (study, chr) that contains a recovered gene, write a checkpoint
# holding the current output rows EXCEPT those genes.
log "seeding checkpoints"
python3 - "$REC" <<'PY' 2>&1 | tee -a "$LOG"
import csv, os, sys, collections
FM="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
RR=os.path.join(FM,"results/susie_coloc_rerun")
bychr=collections.defaultdict(set)
for line in open(sys.argv[1]):
    d,g=line.rstrip("\n").split("\t"); bychr[d.replace("chr","")].add(g)
seeded=tasks=0
for study in sorted(os.listdir(RR)):
    sdir=os.path.join(RR,study)
    if not os.path.isdir(sdir): continue
    for c,genes in bychr.items():
        out=os.path.join(sdir,f"susie_coloc_chr{c}.csv")
        if not os.path.exists(out): continue
        with open(out) as fh:
            rd=csv.DictReader(fh); hdr=rd.fieldnames
            keep=[r for r in rd if r["ensembl"] not in genes]
            total=len(keep)
        # only seed when this file actually contains one of the recovered genes
        with open(out) as fh:
            n_all=sum(1 for _ in csv.DictReader(fh))
        if n_all==total:
            continue
        ck=os.path.join(sdir,f"checkpoint_chr{c}.csv")
        with open(ck,"w",newline="") as fh:
            w=csv.DictWriter(fh,fieldnames=hdr); w.writeheader(); w.writerows(keep)
        seeded+=1; tasks+= (n_all-total)
print(f"  seeded {seeded} checkpoints; {tasks} gene-study pairs will recompute")
PY

log "MERGE + SEED COMPLETE — ready for the targeted COLOC rerun"
