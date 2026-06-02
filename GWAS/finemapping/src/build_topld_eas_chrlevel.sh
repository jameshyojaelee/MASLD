#!/bin/bash -l
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --partition=cpu
#SBATCH --array=1-22
#SBATCH --job-name=build_topld_eas_chr
#SBATCH --output=GWAS/finemapping/logs/build_topld_eas_chr_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/build_topld_eas_chr_%A_%a.err

# Build chr-level PLINK1.9 .bed/.bim/.fam for topld_eas, parallel to topld_eur.
# SuSiEX EAS arm reads chr-level files (10_run_susiex.py:get_eas_ref_prefix lines 352-356),
# NOT per-block. This is required to run SuSiEX with EAS arm = topld_eas instead of 1kg_eas.
#
# Strategy (mirrors backfill_topld_eur_bed.sh, but chr-level instead of per-block):
#   1. Collect union of variant positions across all topld_eas/chr<N>/<bs>.<be>/*.bim files
#   2. Step 1 PLINK extract from 1kg_eas/chr<N>_eas by position range
#   3. Step 2 PLINK extract by SNP-ID matched to topld positions
#   4. Restore topld variant set via concat of per-block .bim files (topld variant IDs/order)

set -o pipefail
module load PLINK/1.9

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
cd "${BASE}"
CHR=${SLURM_ARRAY_TASK_ID}

TOPLD_DIR=data/ld_ref/topld_eas
KG_DIR=data/ld_ref/1kg_eas
KG_PLINK="${KG_DIR}/chr${CHR}_eas"
OUT_PREFIX="${TOPLD_DIR}/chr${CHR}_eas"

echo "[build_topld_eas_chr] chr${CHR} start $(date)"

# Sanity: 1kg_eas chr file must exist
if [ ! -f "${KG_PLINK}.bed" ]; then
  echo "  ERROR: ${KG_PLINK}.bed missing — cannot extract"
  exit 1
fi

# Sanity: topld_eas chr<N>/ subdirectory of blocks must exist
if [ ! -d "${TOPLD_DIR}/chr${CHR}" ]; then
  echo "  ERROR: ${TOPLD_DIR}/chr${CHR} missing — no per-block .bim to anchor against"
  exit 1
fi

# Skip if chr-level PLINK already built (idempotent)
if [ -f "${OUT_PREFIX}.bed" ] && [ -f "${OUT_PREFIX}.bim" ] && [ -f "${OUT_PREFIX}.fam" ]; then
  echo "  SKIP: ${OUT_PREFIX}.{bed,bim,fam} already present"
  exit 0
fi

# Step 0: concat all per-block .bim files for this chr into a single topld bim
# This becomes the "anchor" — defines the variant set + order we want to preserve.
TMP_BIM_CONCAT="${OUT_PREFIX}.topld_concat.bim"
> "${TMP_BIM_CONCAT}"
n_blocks=$(ls -d "${TOPLD_DIR}/chr${CHR}"/*/ 2>/dev/null | wc -l)
for blockdir in "${TOPLD_DIR}/chr${CHR}"/*/; do
  bs_be=$(basename "${blockdir}")
  block_bim="${blockdir}${bs_be}.bim"
  if [ -f "${block_bim}" ]; then
    cat "${block_bim}" >> "${TMP_BIM_CONCAT}"
  fi
done

# Deduplicate by position (column 4); blocks share boundary variants.
# Keep first occurrence — preserves block-order traversal.
awk '!seen[$4]++' "${TMP_BIM_CONCAT}" > "${TMP_BIM_CONCAT}.dedup" && mv "${TMP_BIM_CONCAT}.dedup" "${TMP_BIM_CONCAT}"
n_topld_var=$(wc -l < "${TMP_BIM_CONCAT}")
echo "  topld_eas chr${CHR}: ${n_blocks} blocks, ${n_topld_var} unique-position variants"

# Step 1: range-extract from 1kg_eas (4-col bed: chr, start-1, start, name)
RANGE_BED="${OUT_PREFIX}.range.bed"
awk '{print $1"\t"$4-1"\t"$4"\tv"$4}' "${TMP_BIM_CONCAT}" > "${RANGE_BED}"

STEP1="${OUT_PREFIX}_step1"
if ! plink --bfile "${KG_PLINK}" \
     --extract range "${RANGE_BED}" \
     --make-bed \
     --out "${STEP1}" \
     --allow-extra-chr \
     2>/dev/null; then
  echo "  ERROR: PLINK step 1 failed for chr${CHR}"
  rm -f "${TMP_BIM_CONCAT}" "${RANGE_BED}" "${STEP1}"*
  exit 1
fi
n_step1=$(wc -l < "${STEP1}.bim")
echo "  step1: ${n_step1} variants extracted by position range from 1kg_eas"

# Step 2: identify 1kg_eas SNP IDs whose pos matches topld positions; sub-extract.
# This filters out 1kg_eas positions that fall inside the range but don't match topld
# (e.g., 1kg has more variants than topld in the same window).
KEEP_IDS="${OUT_PREFIX}.keep_ids.txt"
awk 'NR==FNR { pos[$4]=1; next } ($4 in pos) { print $2 }' \
    "${TMP_BIM_CONCAT}" "${STEP1}.bim" > "${KEEP_IDS}"
n_keep=$(wc -l < "${KEEP_IDS}")
echo "  step2: ${n_keep} SNPs match topld positions"

if ! plink --bfile "${STEP1}" \
     --extract "${KEEP_IDS}" \
     --make-bed \
     --out "${OUT_PREFIX}" \
     --allow-extra-chr \
     2>/dev/null; then
  echo "  ERROR: PLINK step 2 failed for chr${CHR}"
  rm -f "${TMP_BIM_CONCAT}" "${RANGE_BED}" "${KEEP_IDS}" "${STEP1}"* "${OUT_PREFIX}.bed" "${OUT_PREFIX}.bim" "${OUT_PREFIX}.fam"
  exit 1
fi

# Final: replace .bim with topld variant IDs (matching per-block .ld matrix ordering).
# Use only the variants that survived the PLINK extract (some 1kg_eas positions may not
# exist; we keep the topld bim subset that intersects).
awk 'NR==FNR { keep[$4]=1; next } ($4 in keep)' \
    "${OUT_PREFIX}.bim" "${TMP_BIM_CONCAT}" > "${OUT_PREFIX}.bim.topld_final"
mv "${OUT_PREFIX}.bim" "${OUT_PREFIX}.bim.kg_eas"
mv "${OUT_PREFIX}.bim.topld_final" "${OUT_PREFIX}.bim"

n_final=$(wc -l < "${OUT_PREFIX}.bim")
echo "  final: ${n_final} variants in ${OUT_PREFIX}.{bed,bim,fam}"

# Cleanup
rm -f "${TMP_BIM_CONCAT}" "${RANGE_BED}" "${KEEP_IDS}" "${STEP1}".{bed,bim,fam,log,nosex} "${OUT_PREFIX}.bim.kg_eas"

echo "[build_topld_eas_chr] chr${CHR} done $(date) — ${n_final}/${n_topld_var} topld variants represented"
