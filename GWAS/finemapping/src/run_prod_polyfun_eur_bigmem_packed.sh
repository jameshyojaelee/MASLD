#!/bin/bash -l
#SBATCH --cpus-per-task=28
#SBATCH --mem=1024G
#SBATCH --time=90:00:00
#SBATCH --partition=bigmem
#SBATCH --job-name=polyfun-susie-coloc
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_polyfun_packed_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_polyfun_packed_%A_%a.err
#
# Phase 9d — packed bigmem variant of run_prod_polyfun_eur.sh
# ============================================================
# Runs TASKS_PER_BATCH parallel SuSiE-COLOC tasks per bigmem job, packing
# the 2 TB / 72-CPU bigmem nodes that would otherwise be idle if the regular
# array allocated one task per node.
#
# Memory accounting (per task: peak ~14 GB):
#   15 × 16 GB = 240 GB allocated  (2 TB node has plenty of headroom)
# CPU accounting (per task: 2 CPUs is enough for data.table; R is mostly serial):
#   15 × 2 = 30 CPUs    (72 CPU node has plenty of headroom)
# Wall time:
#   max(individual COLOC) ≈ 60 min × 1 wait cycle = ~60 min per batch
#   --time=6:00:00 leaves 6× safety margin
#
# To submit (after build dependencies cleared):
#   sbatch --array=0-24 \
#          --dependency=afterany:<build_jobs> \
#          run_prod_polyfun_eur_bigmem_packed.sh
#
# 25 batches × 15 tasks/batch = 375 task slots ≥ 374 actual tasks.
# Idempotent skip: if regular array (job 15596563) writes an output first,
# the packed task sees the file and returns 0 without re-running. So safe to
# run concurrently with the regular array.

set -o pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${BASE}"

GWAS_LIST=(
  2019_31311600_NAFLD_EUR  2020_32298765_NAFLD_EUR  2021_34128465_PDFF_EUR
  2021_34841290_NAFLD_EUR  2021_34957434_PDFF_EUR  2022_36402844_PDFF_EUR
  2023_36280732_NAFLD_deCode_EUR  2023_36280732_NAFLD_Intermountain_EUR
  2023_36280732_NAFLD_UKBB_EUR  UKBB_ALT  UKBB_AST  UKBB_GGT
  FinnGen_NAFLD  FinnGen_NASH  FinnGen_HCC  Ghouse_Cirrhosis  Ghouse_HCC
)

# --- Parameters ---------------------------------------------------------
TASKS_PER_BATCH=${TASKS_PER_BATCH:-7}
TOTAL_TASKS=374
BATCH_ID=${SLURM_ARRAY_TASK_ID:-0}
TASK_START=$(( BATCH_ID * TASKS_PER_BATCH ))
TASK_END=$(( TASK_START + TASKS_PER_BATCH - 1 ))
[ "$TASK_END" -ge "$TOTAL_TASKS" ] && TASK_END=$(( TOTAL_TASKS - 1 ))

if [ "$TASK_START" -ge "$TOTAL_TASKS" ]; then
  echo "[coloc_packed] BATCH_ID=${BATCH_ID}: range ${TASK_START}-${TASK_END} out of bounds, exiting"
  exit 0
fi

echo "============================================================"
echo "[coloc_packed] BATCH_ID=${BATCH_ID}"
echo "[coloc_packed] processing tasks ${TASK_START}..${TASK_END}"
echo "[coloc_packed] start $(date)"
echo "[coloc_packed] node $(hostname)  cpus=${SLURM_CPUS_ON_NODE:-?}  mem=${SLURM_MEM_PER_NODE:-?}M"
echo "============================================================"

export LD_PANEL=polyfun
export COLOC_OUT_SUFFIX="_polyfun"
unset UKBB_LD_DIR EAS_LD_DIR AFR_LD_DIR SAS_LD_DIR

# Per-task processing function
process_task() {
  local IDX=$1
  local GWAS_IDX=$((IDX / 22))
  local CHR=$((IDX % 22 + 1))
  local GWAS=${GWAS_LIST[$GWAS_IDX]}
  local OUT_CSV="${BASE}/GWAS/finemapping/results/susie_coloc_polyfun/${GWAS}/susie_coloc_chr${CHR}.csv"

  # Idempotency: skip if output already exists (regular array may have done it)
  if [ -f "${OUT_CSV}" ]; then
    local SIZE=$(stat -c %s "${OUT_CSV}")
    if [ "${SIZE}" -gt 1000 ]; then
      echo "[task ${IDX} ${GWAS} chr${CHR}] skip — output exists (${SIZE} bytes)"
      return 0
    fi
  fi

  echo "[task ${IDX} ${GWAS} chr${CHR}] start $(date +%H:%M:%S)"
  Rscript GWAS/finemapping/src/06_susie_coloc.R "${GWAS}" "${CHR}" 2>&1 \
    | sed "s/^/[task ${IDX} ${GWAS} chr${CHR}] /"
  local rc=$?
  echo "[task ${IDX} ${GWAS} chr${CHR}] done $(date +%H:%M:%S) rc=${rc}"
  return $rc
}

# Fork all tasks in parallel
pids=()
indices=()
for IDX in $(seq ${TASK_START} ${TASK_END}); do
  echo "  → forking task ${IDX}"
  ( process_task ${IDX} ) &
  pids+=($!)
  indices+=(${IDX})
  # Brief stagger to avoid simultaneous LD-file open thundering herd
  sleep 1
done

# Collect results
n_ok=0
n_skip=0
n_fail=0
for i in "${!pids[@]}"; do
  pid=${pids[$i]}
  idx=${indices[$i]}
  if wait "$pid"; then
    n_ok=$((n_ok+1))
  else
    n_fail=$((n_fail+1))
    echo "[coloc_packed] task ${idx} FAILED (pid ${pid})"
  fi
done

echo "============================================================"
echo "[coloc_packed] BATCH_ID=${BATCH_ID} done"
echo "[coloc_packed] OK=${n_ok}  FAIL=${n_fail}  (of $((TASK_END - TASK_START + 1)))"
echo "[coloc_packed] end $(date)"
echo "============================================================"
