#!/bin/bash
# run_atlas_rebuild_chain.sh
# Full atlas rebuild + honest convergence re-score after the 2026-06-01 modality
# expansion (7 new layers merged in 27a). Linear afterok chain:
#   27a (assemble + merge 7 new layers) -> 75 (causal) -> 217 (stratified)
#   -> 46d (convergence heuristic, unchanged axes) -> 46e (Fisher, +3 new axes)
#   -> 27b (benchmark presets) -> portal JSON
# Submits and returns immediately; monitor with: squeue --name=rebuild -u $USER
set -euo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
LOGS=$BASE/RNA-seq/logs/rebuild; mkdir -p "$LOGS"
cd "$BASE"

sub() {  # name dep cpus mem cmd...
  local name=$1 dep=$2 cpus=$3 mem=$4; shift 4
  local depflag=""; [ "$dep" != "none" ] && depflag="--dependency=afterok:$dep"
  sbatch --parsable --job-name=rebuild --partition=cpu --qos=nslab \
    --cpus-per-task="$cpus" --mem="$mem" --time=48:00:00 \
    --output="$LOGS/${name}_%j.out" --error="$LOGS/${name}_%j.err" \
    $depflag --wrap="cd $BASE && $*"
}

J27a=$(sub 27a   "${START_DEP:-none}" 4 64G "$RS RNA-seq/27a_assemble_evidence_atlas.R")
J75=$( sub 75   "$J27a" 4 48G "$RS RNA-seq/75_integrate_causal_overhaul.R")
J217=$(sub 217  "$J75"  4 48G "$RS RNA-seq/217_stratified_causal_atlas.R")
J46d=$(sub 46d  "$J217" 16 64G "$RS RNA-seq/46d_convergence_evidence.R")
J46e=$(sub 46e  "$J46d" 4 48G "$RS RNA-seq/46e_fisher_combined_test.R")
J27b=$(sub 27b  "$J46e" 4 32G "$RS RNA-seq/27b_benchmark_presets.R")
Jport=$(sub portal "$J27b" 2 16G "$PY scripts/portal/generate_convergence_evidence_json.py")

echo "REBUILD CHAIN SUBMITTED (job-name=rebuild):"
echo "  27a=$J27a -> 75=$J75 -> 217=$J217 -> 46d=$J46d -> 46e=$J46e -> 27b=$J27b -> portal=$Jport"
echo "$J27a $J75 $J217 $J46d $J46e $J27b $Jport" > "$LOGS/last_chain_jobids.txt"
