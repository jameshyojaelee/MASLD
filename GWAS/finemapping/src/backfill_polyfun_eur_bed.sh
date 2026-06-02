#!/bin/bash -l
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --partition=cpu,io
#SBATCH --array=1-22
#SBATCH --job-name=backfill_polyfun_eur_bed
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/backfill_polyfun_eur_bed_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/backfill_polyfun_eur_bed_%A_%a.err

# Phase 8f prep — backfill .bed/.fam per polyfun_eur per-block subdir.
# SuSiEX EUR arm needs per-block .bed/.bim/.fam (10_run_susiex.py:415).
# PolyFun ships only LD matrices (.npz/.gz) — no genotypes — so we extract
# matching variants from 1kg_eur into per-block .bed/.fam.
#
# This version bakes in the lesson from topld_eur: when 1kg has duplicate
# positions (multi-allelic), PLINK extract drops them, leaving .bed shorter
# than the polyfun .bim. We KEEP plink's output bim (don't restore polyfun's)
# AND subset polyfun's .ld matrix to match (using strand-aware allele
# matching with first-match-wins on duplicates — same logic as
# fix_topld_eur_block_size.sh).

set -o pipefail
module load PLINK/1.9 2>/dev/null || true

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
cd "${BASE}"
CHR=${SLURM_ARRAY_TASK_ID}

POLYFUN_DIR=data/ld_ref/polyfun_eur
KG_DIR=data/ld_ref/1kg_eur
KG_PLINK="${KG_DIR}/chr${CHR}_eur"

echo "[backfill_polyfun] chr${CHR} start $(date)"
N_BLOCKS=0; N_OK=0; N_SKIP=0; N_FAIL=0

if [ ! -d "${POLYFUN_DIR}/chr${CHR}" ]; then
  echo "  ERROR: ${POLYFUN_DIR}/chr${CHR} missing — run 8b first"; exit 1
fi

for blockdir in "${POLYFUN_DIR}/chr${CHR}"/*/; do
  bs_be=$(basename "${blockdir}")
  PREFIX="${blockdir}${bs_be}"
  N_BLOCKS=$((N_BLOCKS+1))

  # Skip if .bed/.fam already present and consistent
  if [ -f "${PREFIX}.bed" ] && [ -f "${PREFIX}.fam" ]; then
    NSAM=$(wc -l < "${PREFIX}.fam")
    NVAR=$(wc -l < "${PREFIX}.bim")
    EXP_BED=$(( 3 + NVAR * (NSAM + 3) / 4 ))
    ACT_BED=$(stat -c %s "${PREFIX}.bed")
    if [ "${EXP_BED}" = "${ACT_BED}" ]; then
      N_SKIP=$((N_SKIP+1))
      continue
    fi
  fi

  # Need polyfun bim to anchor variant set
  if [ ! -f "${PREFIX}.bim" ]; then
    echo "  WARNING: no polyfun .bim for ${bs_be} — skipping (8b incomplete?)"
    N_FAIL=$((N_FAIL+1))
    continue
  fi
  if [ ! -f "${PREFIX}.ld" ]; then
    echo "  WARNING: no polyfun .ld for ${bs_be} — skipping (8b incomplete?)"
    N_FAIL=$((N_FAIL+1))
    continue
  fi

  WORK=$(mktemp -d)

  # Step 1: range extract from 1kg_eur
  awk '{print $1"\t"$4-1"\t"$4"\tv"$4}' "${PREFIX}.bim" > "${WORK}/range.bed"
  if ! plink --bfile "${KG_PLINK}" \
       --extract range "${WORK}/range.bed" \
       --make-bed --out "${WORK}/step1" \
       --allow-extra-chr 2>/dev/null; then
    rm -rf "${WORK}"
    N_FAIL=$((N_FAIL+1))
    continue
  fi

  # Step 2: keep 1kg SNPs whose position is in polyfun
  awk 'NR==FNR { pos[$4]=1; next } ($4 in pos) { print $2 }' \
      "${PREFIX}.bim" "${WORK}/step1.bim" > "${WORK}/keep_ids.txt"

  if ! plink --bfile "${WORK}/step1" \
       --extract "${WORK}/keep_ids.txt" \
       --make-bed --out "${WORK}/final" \
       --allow-extra-chr 2>/dev/null || [ ! -f "${WORK}/final.bed" ]; then
    rm -rf "${WORK}"
    N_FAIL=$((N_FAIL+1))
    continue
  fi

  NEW_NVAR=$(wc -l < "${WORK}/final.bim")
  POLYFUN_NVAR=$(wc -l < "${PREFIX}.bim")

  # Build keep_mask using strand-aware allele match + first-match-wins (same as
  # fix_topld_eur_block_size.sh)
  awk '{count[$4]++} END {for (p in count) print p, count[p]}' \
      "${PREFIX}.bim" > "${WORK}/polyfun_pos_counts.txt"
  awk '{
    a=$5; b=$6
    if (a > b) { t=a; a=b; b=t }
    print $4, a, b
  }' "${WORK}/final.bim" > "${WORK}/final_pos_alleles.txt"

  awk '
  function comp(b) {
    if (b == "A") return "T"; if (b == "T") return "A"
    if (b == "C") return "G"; if (b == "G") return "C"
    return b
  }
  function sort2(a, b,    t) {
    if (a > b) { t=a; a=b; b=t }
    return a "/" b
  }
  FILENAME ~ /polyfun_pos_counts/ { polyfun_count[$1] = $2; next }
  FILENAME ~ /final_pos_alleles/ { final_alleles[$1] = sort2($2, $3); next }
  {
    pos = $4; a = $5; b = $6
    polyfun_key = sort2(a, b)
    polyfun_comp_key = sort2(comp(a), comp(b))
    if (!(pos in final_alleles)) { print "0"; next }
    if (polyfun_count[pos] == 1) { print "1"; next }
    if (claimed[pos]) { print "0"; next }
    if (polyfun_key == final_alleles[pos] || polyfun_comp_key == final_alleles[pos]) {
      claimed[pos] = 1; print "1"
    } else { print "0" }
  }
  ' "${WORK}/polyfun_pos_counts.txt" "${WORK}/final_pos_alleles.txt" "${PREFIX}.bim" \
    > "${WORK}/keep_mask.txt"

  N_KEEP=$(grep -c "^1" "${WORK}/keep_mask.txt" || true)
  if [ "${N_KEEP}" != "${NEW_NVAR}" ]; then
    echo "  WARNING: ${bs_be} keep_mask=${N_KEEP} != PLINK extract=${NEW_NVAR}; skipping"
    rm -rf "${WORK}"
    N_FAIL=$((N_FAIL+1))
    continue
  fi

  # Subset polyfun .ld using keep_mask
  awk -v maskfile="${WORK}/keep_mask.txt" '
  BEGIN {
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
  }' "${PREFIX}.ld" > "${WORK}/new.ld"

  # Backup and replace
  cp "${PREFIX}.bim" "${PREFIX}.bim.bak_pre_backfill" 2>/dev/null
  cp "${PREFIX}.ld"  "${PREFIX}.ld.bak_pre_backfill"  2>/dev/null
  cp "${WORK}/final.bim" "${PREFIX}.bim"
  cp "${WORK}/final.bed" "${PREFIX}.bed"
  cp "${WORK}/final.fam" "${PREFIX}.fam"
  cp "${WORK}/new.ld"    "${PREFIX}.ld"

  # Verify final consistency
  FINAL_NVAR=$(wc -l < "${PREFIX}.bim")
  FINAL_LD_R=$(wc -l < "${PREFIX}.ld")
  FINAL_LD_C=$(awk 'NR==1{print NF; exit}' "${PREFIX}.ld")
  FINAL_BED=$(stat -c %s "${PREFIX}.bed")
  FINAL_FAM=$(wc -l < "${PREFIX}.fam")
  EXP_BED=$(( 3 + FINAL_NVAR * (FINAL_FAM + 3) / 4 ))

  if [ "${FINAL_NVAR}" = "${FINAL_LD_R}" ] && \
     [ "${FINAL_LD_R}" = "${FINAL_LD_C}" ] && \
     [ "${EXP_BED}" = "${FINAL_BED}" ]; then
    rm -f "${PREFIX}.bim.bak_pre_backfill" "${PREFIX}.ld.bak_pre_backfill"
    N_OK=$((N_OK+1))
  else
    cp "${PREFIX}.bim.bak_pre_backfill" "${PREFIX}.bim" 2>/dev/null
    cp "${PREFIX}.ld.bak_pre_backfill"  "${PREFIX}.ld"  2>/dev/null
    rm -f "${PREFIX}.bed" "${PREFIX}.fam" "${PREFIX}.bim.bak_pre_backfill" "${PREFIX}.ld.bak_pre_backfill"
    N_FAIL=$((N_FAIL+1))
  fi
  rm -rf "${WORK}"
done

echo "[backfill_polyfun] chr${CHR} done: ${N_OK} built, ${N_SKIP} skipped, ${N_FAIL} failed ($(date))"
