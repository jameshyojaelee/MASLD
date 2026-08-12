#!/bin/bash
# run_atlas_rebuild_chain.sh
# Full atlas rebuild + honest convergence re-score after the 2026-06-01 modality
# expansion (7 new layers merged in 27a). Linear afterok chain:
#   27a (assemble + merge 7 new layers) -> 75 (causal) -> 217 (stratified)
#   -> 46d (convergence heuristic, unchanged axes) -> 46e (Fisher, +3 new axes)
#   -> 27b (benchmark presets). The historical rank-first portal export is
#   disabled by default and is not the MASLD Gene Catalog release.
# Submits and returns immediately; monitor with: squeue --name=rebuild -u $USER
set -euo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
LOGS=$BASE/RNA-seq/logs/rebuild; mkdir -p "$LOGS"
cd "$BASE"

sub() {  # name dep cpus mem partition qos cmd...
  local name=$1 dep=$2 cpus=$3 mem=$4 part=$5 qos=$6; shift 6
  local depflag=""; [ "$dep" != "none" ] && depflag="--dependency=afterok:$dep"
  sbatch --parsable --job-name=rebuild --partition="$part" --qos="$qos" \
    --cpus-per-task="$cpus" --mem="$mem" --time=48:00:00 \
    --output="$LOGS/${name}_%j.out" --error="$LOGS/${name}_%j.err" \
    $depflag --wrap="cd $BASE && $*"
}

# C2 migration (2026-06-15): added 45a (recomputes sources_active from the migrated
# bulk_* human-bulk channel — required after the dream->bulk swap); 46d on bigmem
# (perm10k FDR on the 27,187-gene C2 universe needs the RAM). START_DEP lets a
# verified standalone 27a feed the chain from 45a (set START_DEP=<27a_jobid> and
# comment out the J27a line).
J27a=$(sub 27a   "${START_DEP:-none}" 4 64G  cpu    nslab "$RS RNA-seq/27a_assemble_evidence_atlas.R")
J45a=$(sub 45a   "$J27a" 4 48G  cpu    nslab "$RS RNA-seq/45a_integrate_all_sources.R")
J75=$( sub 75    "$J45a" 4 48G  cpu    nslab "$RS RNA-seq/75_integrate_causal_overhaul.R")
J217=$(sub 217   "$J75"  4 48G  cpu    nslab "$RS RNA-seq/217_stratified_causal_atlas.R")
J46d=$(sub 46d   "$J217" 16 192G cpu    nslab "$RS RNA-seq/46d_convergence_evidence.R")  # cpu, not bigmem: bigmem enforces a 500G min (QOSMinMemory); 192G fits the cpu 210G cap
J46e=$(sub 46e   "$J46d" 4 48G  cpu    nslab "$RS RNA-seq/46e_fisher_combined_test.R")
J27b=$(sub 27b   "$J46e" 4 32G  cpu    nslab "$RS RNA-seq/27b_benchmark_presets.R")
Jport="SKIPPED_LEGACY_PORTAL"
if [[ "${ALLOW_LEGACY_PORTAL_REBUILD:-false}" == "true" ]]; then
  Jport=$(sub portal "$J27b" 2 16G cpu nslab \
    "ALLOW_LEGACY_PORTAL_REBUILD=true $PY scripts/portal/generate_convergence_evidence_json.py")
fi

echo "REBUILD CHAIN SUBMITTED (job-name=rebuild):"
echo "  27a=$J27a -> 45a=$J45a -> 75=$J75 -> 217=$J217 -> 46d=$J46d -> 46e=$J46e -> 27b=$J27b; legacy portal=$Jport"
echo "$J27a $J45a $J75 $J217 $J46d $J46e $J27b $Jport" > "$LOGS/last_chain_jobids.txt"
