#!/usr/bin/env bash
# Creates a fresh direct-only preview root (symlinked inputs) and points run_preview_direct.sbatch at it. Prints the root.
set -euo pipefail
P=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design; RES=$P/GWAS/finemapping/results/alphagenome_atlas; MAIN=$RES/run-20260909T153939Z
TS=$(date -u +%Y%m%dT%H%M%SZ); PREV=$RES/preview-direct-only-$TS
mkdir -p $PREV/tables $PREV/raw $PREV/logs $PREV/validation
for f in $MAIN/tables/*; do ln -s $f $PREV/tables/; done
for d in atlas_direct pilot scorer_metadata atlas_gate; do ln -s $MAIN/raw/$d $PREV/raw/$d; done
echo "pipeline check on the direct-universe archive only (13,690 variants) + gate archive; NOT a result; final tables come from run-20260909T153939Z" > $PREV/PREVIEW_NOTE.txt
sed -i -E "s#preview-direct-only-[0-9TZ]+#preview-direct-only-$TS#g" $P/scripts/analysis/alphagenome_atlas/run_preview_direct.sbatch
echo $PREV
