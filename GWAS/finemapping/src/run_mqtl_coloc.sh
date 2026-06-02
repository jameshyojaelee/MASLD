#!/bin/bash
#SBATCH --job-name=coloc
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --array=0-19
#SBATCH --output=GWAS/finemapping/logs/60b_mqtl_coloc_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/60b_mqtl_coloc_%A_%a.err
# ---------------------------------------------------------------------------
# run_mqtl_coloc.sh
# SLURM array for the metabolite/lipid-QTL COLOC layer (ABF-only first pass).
#
# Each array task runs BOTH arms (Chen + Ottensmann) for its shard index,
# with N_SHARDS = array size. coloc.abf, MASLD GWAS credible-set loci only.
#
# Prereq: 60a download job done (sentinel data/external/chen2023_mqtl/.download_complete)
#         60_prep_mqtl_loci.R already produced masld_credset_loci.tsv.
#
# Submit:  cd <root> && sbatch GWAS/finemapping/src/run_mqtl_coloc.sh
# Combine: micromamba run -n rnaseq Rscript GWAS/finemapping/src/61b_combine_metabolite_coloc.R
# ---------------------------------------------------------------------------
set -uo pipefail
BASE_DIR="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$BASE_DIR" || exit 1

N_SHARDS=$(( ${SLURM_ARRAY_TASK_COUNT:-20} ))
SHARD_ID=${SLURM_ARRAY_TASK_ID:-0}

echo "=== mQTL COLOC shard $SHARD_ID / $N_SHARDS | $(date) ==="

# Wait (politely) for the download sentinel if the array launched before the
# download finished. Bounded wait so a stuck array task eventually exits.
SENTINEL="$BASE_DIR/data/external/chen2023_mqtl/.download_complete"
waited=0
while [[ ! -f "$SENTINEL" && $waited -lt 86400 ]]; do
  echo "  waiting for download sentinel ($waited s)..."; sleep 120; waited=$((waited+120))
done
if [[ ! -f "$SENTINEL" ]]; then echo "ERROR: download sentinel never appeared"; exit 1; fi

for SRC in chen ottensmann; do
  echo "--- arm: $SRC ---"
  QTL_SOURCE=$SRC SHARD_ID=$SHARD_ID N_SHARDS=$N_SHARDS \
    micromamba run -n rnaseq Rscript GWAS/finemapping/src/60b_mqtl_coloc.R
  echo "--- $SRC arm exit: $? ---"
done

echo "=== shard $SHARD_ID done | $(date) ==="
