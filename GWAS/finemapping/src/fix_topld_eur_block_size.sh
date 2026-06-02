#!/bin/bash
# fix_topld_eur_block_size.sh — fix .bim/.ld/.bed inconsistency in topld_eur blocks.
#
# Symptom: SuSiEX errors with "Invalid .bed file size (expected N bytes)."
# Root cause: backfill_topld_eur_bed.sh restored the full TOP-LD .bim AFTER
# PLINK step 2 extracted only the 1kg-matched subset. The .bed has fewer
# variants than .bim claims (e.g., when 1kg has duplicate-position variants
# at a topld position, plink --extract drops both → .bed missing those).
#
# Fix strategy:
#   1. Re-run the 2-step PLINK extract (range → ID match) WITHOUT restoring
#      the topld bim. This produces a .bim consistent with .bed.
#   2. Subset the .ld matrix to keep only rows/cols whose POSITION is
#      represented in the new .bim. The .ld came from TOP-LD original release;
#      its order matches topld's full variant set, so we drop the rows/cols
#      corresponding to topld variants that didn't survive the 1kg extract.
#
# Usage: bash fix_topld_eur_block_size.sh chr21 44506268.46177105

set -euo pipefail
module load PLINK/1.9 2>/dev/null || true

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
cd "${BASE}"

CHR="$1"
BLK="$2"
TOPLD_DIR="data/ld_ref/topld_eur"
KG_DIR="data/ld_ref/1kg_eur"
KG_PLINK="${KG_DIR}/chr${CHR#chr}_eur"
BD="${TOPLD_DIR}/${CHR}/${BLK}"
PREFIX="${BD}/${BLK}"

if [ ! -f "${PREFIX}.bim" ] || [ ! -f "${PREFIX}.bed" ] || [ ! -f "${PREFIX}.ld" ]; then
  echo "ERROR: missing files in ${BD}"; exit 1
fi

NSAM=$(wc -l < "${PREFIX}.fam")
NVAR_BIM=$(wc -l < "${PREFIX}.bim")
BYTES_PER_VAR=$(( (NSAM + 3) / 4 ))
EXPECTED_BED=$(( 3 + NVAR_BIM * BYTES_PER_VAR ))
ACTUAL_BED=$(stat -c %s "${PREFIX}.bed")

echo "Block: ${CHR}/${BLK}"
echo "  bim variants: ${NVAR_BIM}"
echo "  fam samples : ${NSAM}"
echo "  expected bed: ${EXPECTED_BED} bytes"
echo "  actual bed  : ${ACTUAL_BED} bytes"

if [ "${EXPECTED_BED}" = "${ACTUAL_BED}" ]; then
  echo "  OK — sizes match, nothing to fix"; exit 0
fi

# Save originals
cp "${PREFIX}.bim" "${PREFIX}.bim.bak_pre_fix"
cp "${PREFIX}.bed" "${PREFIX}.bed.bak_pre_fix"
cp "${PREFIX}.fam" "${PREFIX}.fam.bak_pre_fix"
cp "${PREFIX}.ld"  "${PREFIX}.ld.bak_pre_fix"

# Re-run the 2-step PLINK extract from 1kg_eur, this time keeping plink's bim
WORK=$(mktemp -d)
trap "rm -rf ${WORK}" EXIT

# Range bed: chr start-1 start name (4 cols)
awk '{print $1"\t"$4-1"\t"$4"\tv"$4}' "${PREFIX}.bim.bak_pre_fix" > "${WORK}/range.bed"

# Step 1: range extract from 1kg_eur
if ! plink --bfile "${KG_PLINK}" \
     --extract range "${WORK}/range.bed" \
     --make-bed \
     --out "${WORK}/step1" \
     --allow-extra-chr 2>/dev/null; then
  echo "ERROR: PLINK step 1 failed"; exit 1
fi

# Step 2: keep only 1kg variants whose position is in topld bim
awk 'NR==FNR { pos[$4]=1; next } ($4 in pos) { print $2 }' \
    "${PREFIX}.bim.bak_pre_fix" "${WORK}/step1.bim" > "${WORK}/keep_ids.txt"

# Final extract; THIS time we KEEP plink's output bim (no topld restore).
# PLINK 1.9 silently drops duplicates during --extract, which is exactly what
# created the original mismatch — bim claimed N+1 but bed got N.
if ! plink --bfile "${WORK}/step1" \
     --extract "${WORK}/keep_ids.txt" \
     --make-bed \
     --out "${WORK}/final" \
     --allow-extra-chr 2>/dev/null; then
  echo "ERROR: PLINK step 2 failed"; exit 1
fi

NEW_NVAR=$(wc -l < "${WORK}/final.bim")
echo "  PLINK final extract: ${NEW_NVAR} variants (matches .bed)"

# Verify the new .bed size matches new bim
NEW_EXPECTED=$(( 3 + NEW_NVAR * BYTES_PER_VAR ))
NEW_ACTUAL=$(stat -c %s "${WORK}/final.bed")
if [ "${NEW_EXPECTED}" != "${NEW_ACTUAL}" ]; then
  echo "ERROR: PLINK output is itself inconsistent (expected ${NEW_EXPECTED}, got ${NEW_ACTUAL})"
  exit 1
fi

# Build keep_mask using a position-then-allele strategy that handles strand flips.
#
# Logic per topld bim row:
#   - If position NOT in new_bim → drop (1kg lacks this position).
#   - If position appears EXACTLY ONCE in topld bim → keep (matches new_bim
#     regardless of strand convention).
#   - If position appears MULTIPLE TIMES in topld bim (e.g., true biallelic
#     duplicates at the same site, like A/G + T/C strand-related) → keep ONLY
#     the topld row whose sorted-alleles match new_bim's sorted-alleles (or
#     the strand-complement of those alleles); drop the others.
#
# Strand complement: A↔T, C↔G.

# Pass 1: count topld occurrences per position
awk '{count[$4]++} END {for (p in count) print p, count[p]}' \
    "${PREFIX}.bim.bak_pre_fix" > "${WORK}/topld_pos_counts.txt"

# Pass 2: read final bim — per position, capture sorted alleles
awk '{
  a=$5; b=$6
  if (a > b) { t=a; a=b; b=t }
  print $4, a, b
}' "${WORK}/final.bim" > "${WORK}/final_pos_alleles.txt"

# Pass 3: decision per topld row, dispatched on FILENAME.
# When topld has duplicate-position rows AND multiple match the new bim's
# alleles (direct + strand-comp both match), keep only the FIRST direct match;
# remaining matches at same pos get dropped. This breaks the 1-row-per-position
# tie naturally.
awk '
function comp(b) {
  if (b == "A") return "T"
  if (b == "T") return "A"
  if (b == "C") return "G"
  if (b == "G") return "C"
  return b
}
function sort2(a, b,    t) {
  if (a > b) { t=a; a=b; b=t }
  return a "/" b
}
FILENAME ~ /topld_pos_counts/ { topld_count[$1] = $2; next }
FILENAME ~ /final_pos_alleles/ { final_alleles[$1] = sort2($2, $3); next }
{
  pos = $4
  a = $5; b = $6
  topld_key = sort2(a, b)
  topld_comp_key = sort2(comp(a), comp(b))

  if (!(pos in final_alleles)) {
    print "0"
    next
  }
  if (topld_count[pos] == 1) {
    print "1"
    next
  }
  # Duplicate position. First-match-wins: keep only ONE topld row per position.
  if (claimed[pos]) {
    print "0"
    next
  }
  # Prefer direct match over strand-comp match. If neither matches, drop.
  if (topld_key == final_alleles[pos]) {
    claimed[pos] = 1
    print "1"
  } else if (topld_comp_key == final_alleles[pos]) {
    claimed[pos] = 1
    print "1"
  } else {
    print "0"
  }
}
' "${WORK}/topld_pos_counts.txt" "${WORK}/final_pos_alleles.txt" "${PREFIX}.bim.bak_pre_fix" > "${WORK}/keep_mask.txt"

N_KEEP=$(grep -c "^1" "${WORK}/keep_mask.txt" || true)
N_DROP=$(grep -c "^0" "${WORK}/keep_mask.txt" || true)
echo "  keep_mask: ${N_KEEP} keep, ${N_DROP} drop (sums to ${NVAR_BIM})"

if [ "${N_KEEP}" != "${NEW_NVAR}" ]; then
  echo "ERROR: keep_mask count (${N_KEEP}) != new bim count (${NEW_NVAR})"
  exit 1
fi

# Subset .ld using keep_mask
awk -v maskfile="${WORK}/keep_mask.txt" '
BEGIN {
  # Read mask
  i=0
  while ((getline line < maskfile) > 0) { i++; mask[i] = line }
  close(maskfile)
}
{
  if (mask[NR] != "1") next
  out=""
  for (j=1; j<=NF; j++) {
    if (mask[j] != "1") continue
    out = (out=="") ? $j : out " " $j
  }
  print out
}' "${PREFIX}.ld.bak_pre_fix" > "${WORK}/new.ld"

# Replace files
cp "${WORK}/final.bim" "${PREFIX}.bim"
cp "${WORK}/final.bed" "${PREFIX}.bed"
cp "${WORK}/final.fam" "${PREFIX}.fam"
cp "${WORK}/new.ld"    "${PREFIX}.ld"

# Verify final state
FINAL_NVAR_BIM=$(wc -l < "${PREFIX}.bim")
FINAL_LD_ROWS=$(wc -l < "${PREFIX}.ld")
FINAL_LD_COLS=$(awk 'NR==1{print NF; exit}' "${PREFIX}.ld")
FINAL_BED_SIZE=$(stat -c %s "${PREFIX}.bed")
FINAL_EXPECTED=$(( 3 + FINAL_NVAR_BIM * BYTES_PER_VAR ))

echo "  AFTER fix:"
echo "    bim: ${FINAL_NVAR_BIM} rows"
echo "    ld : ${FINAL_LD_ROWS} rows × ${FINAL_LD_COLS} cols"
echo "    bed: ${FINAL_BED_SIZE} bytes (expected ${FINAL_EXPECTED})"

if [ "${FINAL_NVAR_BIM}" = "${FINAL_LD_ROWS}" ] && \
   [ "${FINAL_LD_ROWS}" = "${FINAL_LD_COLS}" ] && \
   [ "${FINAL_EXPECTED}" = "${FINAL_BED_SIZE}" ]; then
  echo "  ✓ FIXED"
else
  echo "  ✗ INCONSISTENT — restoring backups"
  cp "${PREFIX}.bim.bak_pre_fix" "${PREFIX}.bim"
  cp "${PREFIX}.bed.bak_pre_fix" "${PREFIX}.bed"
  cp "${PREFIX}.fam.bak_pre_fix" "${PREFIX}.fam"
  cp "${PREFIX}.ld.bak_pre_fix"  "${PREFIX}.ld"
  exit 1
fi
