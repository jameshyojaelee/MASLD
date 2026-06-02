#!/bin/bash -l
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --partition=cpu
#SBATCH --array=1-22
#SBATCH --job-name=backfill_topld_eur_bed
#SBATCH --output=GWAS/finemapping/logs/backfill_topld_eur_bed_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/backfill_topld_eur_bed_%A_%a.err

# Backfill .bed / .fam per topld_eur block. build_topld_blocks.R produced .ld
# and .bim only (variant IDs are TOP-LD's rsIDs). SuSiEX needs .bed/.fam.
# Strategy: extract by position range from 1kg_eur per-block PLINK file (which
# has .bed/.bim/.fam). After PLINK extract, restore topld's original .bim so
# rsIDs are preserved (the .ld matrix order is keyed by topld's variant set).

set -o pipefail
module load PLINK/1.9

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
cd "${BASE}"
CHR=${SLURM_ARRAY_TASK_ID}

TOPLD_DIR=data/ld_ref/topld_eur
KG_DIR=data/ld_ref/1kg_eur
KG_PLINK="${KG_DIR}/chr${CHR}_eur"
BLOCKS_FILE="${TOPLD_DIR}/approx_LD_blocks.txt"

echo "[backfill_topld] chr${CHR} start $(date)"
N_BLOCKS=0; N_OK=0; N_SKIP=0; N_FAIL=0

while read -r bchr bs be; do
  [ "${bchr}" = "chr" ] && continue
  [ "${bchr}" != "${CHR}" ] && continue
  N_BLOCKS=$((N_BLOCKS+1))
  BLOCK_DIR="${TOPLD_DIR}/chr${CHR}/${bs}.${be}"
  BLOCK_PREFIX="${BLOCK_DIR}/${bs}.${be}"

  # Skip if .bed/.fam already present
  if [ -f "${BLOCK_PREFIX}.bed" ] && [ -f "${BLOCK_PREFIX}.fam" ]; then
    N_SKIP=$((N_SKIP+1))
    continue
  fi

  # Need topld bim (positional reference)
  if [ ! -f "${BLOCK_PREFIX}.bim" ]; then
    echo "  WARNING: no topld .bim for ${bs}.${be} — skipping"
    N_FAIL=$((N_FAIL+1))
    continue
  fi

  # Save the topld bim aside
  cp "${BLOCK_PREFIX}.bim" "${BLOCK_PREFIX}.bim.topld_orig"

  # Build position-range file (chr, start-1, start, name) — PLINK requires 4 columns
  awk '{print $1"\t"$4-1"\t"$4"\tv"$4}' "${BLOCK_PREFIX}.bim.topld_orig" > "${BLOCK_PREFIX}.range.bed"

  # Step 1: plink --extract range to get all 1kg_eur variants in topld block window
  STEP1="${BLOCK_PREFIX}_step1"
  if ! plink --bfile "${KG_PLINK}" \
       --extract range "${BLOCK_PREFIX}.range.bed" \
       --make-bed \
       --out "${STEP1}" \
       --allow-extra-chr \
       2>/dev/null; then
    mv "${BLOCK_PREFIX}.bim.topld_orig" "${BLOCK_PREFIX}.bim" 2>/dev/null
    rm -f "${BLOCK_PREFIX}.range.bed" "${STEP1}"*
    N_FAIL=$((N_FAIL+1))
    continue
  fi

  # Step 2: identify 1kg_eur SNP IDs that match topld positions, then sub-extract.
  # build_topld_blocks.R uses half-open [start, stop) and may drop a few variants
  # that 1kg_eur keeps (different boundary handling, lifted-over positions, etc.).
  # This step ensures the final .bed/.bim variant set EXACTLY matches topld's bim
  # so it stays consistent with the precomputed .ld matrix.
  awk 'NR==FNR { pos[$4]=1; next } ($4 in pos) { print $2 }' \
      "${BLOCK_PREFIX}.bim.topld_orig" "${STEP1}.bim" > "${BLOCK_PREFIX}.keep_ids.txt"

  if plink --bfile "${STEP1}" \
      --extract "${BLOCK_PREFIX}.keep_ids.txt" \
      --make-bed \
      --out "${BLOCK_PREFIX}" \
      --allow-extra-chr \
      2>/dev/null && [ -f "${BLOCK_PREFIX}.bed" ]; then

    # Restore topld bim (rsIDs / topld variant set, matches .ld order)
    mv "${BLOCK_PREFIX}.bim.topld_orig" "${BLOCK_PREFIX}.bim"
    rm -f "${BLOCK_PREFIX}.range.bed" "${BLOCK_PREFIX}.keep_ids.txt" "${STEP1}"*
    N_OK=$((N_OK+1))
  else
    mv "${BLOCK_PREFIX}.bim.topld_orig" "${BLOCK_PREFIX}.bim" 2>/dev/null
    rm -f "${BLOCK_PREFIX}.range.bed" "${BLOCK_PREFIX}.keep_ids.txt" "${STEP1}"*
    N_FAIL=$((N_FAIL+1))
  fi
done < "${BLOCKS_FILE}"

echo "[backfill_topld] chr${CHR} done: ${N_OK} built, ${N_SKIP} pre-existing, ${N_FAIL} failed ($(date))"
