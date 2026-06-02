#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --partition=cpu
#SBATCH --array=0-7
#SBATCH --job-name=polyfun_dl
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/polyfun_dl_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/polyfun_dl_%A_%a.err

# Bulk download of PolyFun precomputed UKBB LD bucket
# (broad-alkesgroup-ukbb-ld/UKBB_LD/), 5540 paired files (.npz + .gz), ~2.85 TB.
#
# Strategy:
#   - Pre-chunked into 8 chunks (~693 keys each) at chunks/chunk_NN
#   - Each array task processes one chunk via xargs -P 4 wget -c
#   - --partition=cpu — hundreds of cores available; downloads are NETWORK-bound
#     (S3 throughput), not local-disk-bound, so io partition's optimization
#     doesn't help. cpu has way more headroom than io's 48-core ceiling.
#   - --cpus-per-task=4 keeps memory/CPU footprint small so 8 tasks can run
#     concurrent if any slots are free anywhere on cpu partition.
#   - Default nslab QOS — concurrent-job limit is much higher than interactive QOS.
#
# USAGE:
#   Default (full bucket, 8 chunks, ~30-60 min wall):
#     sbatch /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src/download_polyfun_ld.sh
#
#   Subset (e.g., chunks 4-7 only — when another user already ran 0-3):
#     sbatch --array=4-7 /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src/download_polyfun_ld.sh
#
# SAFETY: Do NOT have two users running the SAME array index concurrently — wget -c
# is not lock-safe across writers and will corrupt partial files. Coordinate via
# --array=N-M so users cover disjoint chunks.
#
# Resume support: wget -c continues partial downloads; safe to re-run after a job dies.
# Output dir + logs dir are group-writable (nslab, mode 775).

set -o pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
cd "${BASE}"

CHUNK_ID=$(printf "%02d" "${SLURM_ARRAY_TASK_ID}")
CHUNK_FILE="data/ld_ref/polyfun_eur/chunks/chunk_${CHUNK_ID}"
DEST_DIR="data/ld_ref/polyfun_eur/raw"
S3_BASE="https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com"

if [ ! -f "${CHUNK_FILE}" ]; then
  echo "ERROR: chunk file ${CHUNK_FILE} missing"
  exit 1
fi

mkdir -p "${DEST_DIR}"

n_files=$(wc -l < "${CHUNK_FILE}")
echo "[polyfun_dl] task ${SLURM_ARRAY_TASK_ID}: ${n_files} files in chunk ${CHUNK_ID}  start $(date)"

# Pre-flight: count how many files already exist (idempotent resume)
n_existing=0
while read -r key; do
  fname=$(basename "${key}")
  if [ -s "${DEST_DIR}/${fname}" ]; then
    n_existing=$((n_existing+1))
  fi
done < "${CHUNK_FILE}"
echo "  ${n_existing}/${n_files} already present (skip via wget -c)"

# Parallel download via xargs -P 8 (8 concurrent wget streams within this task)
# wget -c: continue partial downloads. -q: quiet (suppress per-file progress)
# Output goes to DEST_DIR via -P. -nv: non-verbose summary per file.
awk -v base="${S3_BASE}" '{print base"/"$0}' "${CHUNK_FILE}" \
  | xargs -P 4 -I {} wget -c -nv -t 5 --waitretry=10 -P "${DEST_DIR}" {} 2>&1 \
  | tail -50

# Post-flight verification
n_now=0
n_missing=0
> /tmp/missing_${CHUNK_ID}.txt
while read -r key; do
  fname=$(basename "${key}")
  if [ -s "${DEST_DIR}/${fname}" ]; then
    n_now=$((n_now+1))
  else
    n_missing=$((n_missing+1))
    echo "${key}" >> /tmp/missing_${CHUNK_ID}.txt
  fi
done < "${CHUNK_FILE}"

echo ""
echo "[polyfun_dl] task ${SLURM_ARRAY_TASK_ID}: complete $(date)"
echo "  files OK: ${n_now}/${n_files}"
if [ "${n_missing}" -gt 0 ]; then
  echo "  MISSING: ${n_missing} files (see /tmp/missing_${CHUNK_ID}.txt)"
  head -20 /tmp/missing_${CHUNK_ID}.txt
  exit 2
fi

echo "  chunk_${CHUNK_ID} OK"
