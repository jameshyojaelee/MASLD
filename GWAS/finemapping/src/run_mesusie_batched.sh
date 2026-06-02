#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --mem=280G
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem
#SBATCH --job-name=mesusie_batch
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/mesusie_batch_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/mesusie_batch_%A_%a.err
#
# Batched MESuSiE wrapper — runs 5 loci in parallel within one bigmem task,
# maximizing RAM usage on the 360GB bigmem nodes.
#
# Per-locus peak memory ≈ 48 GB → 5 parallel × 48 GB = 240 GB allocated;
# request 280 GB to leave 40 GB OOM headroom.
#
# Required env vars (set via sbatch --export=ALL,...):
#   MESUSIE_RESULTS_SUFFIX  (e.g., _1kg / _topld / _topld_full)
#   LD_PANEL                (1kg or topld)
#   LOCI_PER_BATCH          (default 5)
#   LOCI_OFFSET             (start offset, e.g., 38 for 1kg-tail; 0 for fresh)
# Optional:
#   EAS_LD_DIR              (pin EAS arm; e.g., 1kg_eas for topld-EUR-only swap)

set -eo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

# Required env
: ${MESUSIE_RESULTS_SUFFIX:?must be set (e.g., _1kg, _topld, _topld_full)}
: ${LD_PANEL:?must be set (1kg or topld)}

LOCI_PER_BATCH=${LOCI_PER_BATCH:-5}
LOCI_OFFSET=${LOCI_OFFSET:-0}
N_TOTAL_LOCI=$(( $(wc -l < results/susiex/shared_loci.csv) - 1 ))
BATCH_ID=${SLURM_ARRAY_TASK_ID:-1}
LOCUS_START=$(( (BATCH_ID - 1) * LOCI_PER_BATCH + LOCI_OFFSET + 1 ))
LOCUS_END=$(( LOCUS_START + LOCI_PER_BATCH - 1 ))
[ "$LOCUS_END" -gt "$N_TOTAL_LOCI" ] && LOCUS_END=$N_TOTAL_LOCI

echo "[mesusie_batch] BATCH_ID=${BATCH_ID}  loci ${LOCUS_START}-${LOCUS_END} of ${N_TOTAL_LOCI}"
echo "[mesusie_batch] suffix=${MESUSIE_RESULTS_SUFFIX}  LD_PANEL=${LD_PANEL}  EAS_LD_DIR=${EAS_LD_DIR:-<unset>}"
echo "[mesusie_batch] start $(date)"

# Launch parallel R workers
pids=()
for LOCUS_ROW in $(seq ${LOCUS_START} ${LOCUS_END}); do
  echo "  → forking LOCUS_ROW=${LOCUS_ROW}"
  (
    LOCUS_ROW=${LOCUS_ROW} Rscript src/10b_run_mesusie.R 2>&1 \
      | sed "s/^/[locus ${LOCUS_ROW}] /"
  ) &
  pids+=($!)
done

# Wait for all workers; collect exit statuses
n_ok=0; n_fail=0
for pid in "${pids[@]}"; do
  if wait "$pid"; then
    n_ok=$((n_ok+1))
  else
    n_fail=$((n_fail+1))
    echo "  WARNING: worker pid ${pid} exited non-zero"
  fi
done

echo "[mesusie_batch] BATCH_ID=${BATCH_ID} done  ${n_ok} OK, ${n_fail} fail  $(date)"
