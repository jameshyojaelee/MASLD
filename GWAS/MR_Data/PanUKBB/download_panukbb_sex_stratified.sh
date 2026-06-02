#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --job-name=panukbb_sex_dl
#SBATCH --output=GWAS/MR_Data/PanUKBB/sex_stratified/logs/download_%j.out
#SBATCH --error=GWAS/MR_Data/PanUKBB/sex_stratified/logs/download_%j.err

# ---------------------------------------------------------------------------
# Download sex-stratified UK Biobank GWAS for liver enzymes (ALT/AST/GGT).
#
# IMPORTANT: Pan-UKBB v0.4 does NOT publish sex-stratified biomarker GWAS
# (only pheno_sex == both_sexes in the manifest). We therefore fall back to
# the Neale Lab Round 2 UKBB GWAS, which DOES publish per-sex files (EUR-only
# but that's what we need against the PolyFun EUR LD reference).
#
# Source bucket: broad-ukb-sumstats-us-east-1
# Path pattern:  round2/additive-tsvs/{phenocode}_irnt.gwas.imputed_v3.{sex}.varorder.tsv.bgz
# Format:        variant (chr:pos:ref:alt hg19), minor_allele, minor_AF, beta, se, pval, n_complete_samples
# ---------------------------------------------------------------------------

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

DEST="GWAS/MR_Data/PanUKBB/sex_stratified"
mkdir -p "${DEST}/logs"

BUCKET="https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs"

declare -A PHENOCODES=( [ALT]=30620 [AST]=30650 [GGT]=30730 )

echo "=== Downloading Neale Lab Round 2 sex-stratified GWAS ==="
echo "Start: $(date)"
echo ""

for trait in ALT AST GGT; do
  pc=${PHENOCODES[$trait]}
  for sex in female male; do
    fn="${pc}_irnt.gwas.imputed_v3.${sex}.varorder.tsv.bgz"
    url="${BUCKET}/${fn}"
    out="${DEST}/${fn}"
    if [ -s "${out}" ]; then
      sz=$(stat -c%s "${out}")
      if [ "${sz}" -gt 400000000 ]; then
        echo "  SKIP (already downloaded): ${fn} (${sz} bytes)"
        continue
      else
        echo "  RE-DOWNLOAD (partial): ${fn}"
        rm -f "${out}"
      fi
    fi
    echo "--- Downloading ${trait} ${sex} (phenocode ${pc}) ---"
    echo "  URL: ${url}"
    curl --fail --silent --show-error --location \
         --output "${out}" \
         "${url}" \
    || { echo "  FAILED: ${fn}"; rm -f "${out}"; exit 2; }
    sz=$(stat -c%s "${out}")
    echo "  Saved: ${out} (${sz} bytes)"
  done
done

# Also fetch the variants annotation file (chr:pos:ref:alt -> rsid + hg19 coords)
VAR_FN="variants.tsv.bgz"
if [ ! -s "${DEST}/${VAR_FN}" ]; then
  echo ""
  echo "--- Downloading variants annotation ---"
  # Round 2 variant manifest has rsid + chr + pos
  curl --fail --silent --show-error --location \
       --output "${DEST}/${VAR_FN}" \
       "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz" \
  || { echo "  WARNING: variants annotation download failed (non-fatal)"; rm -f "${DEST}/${VAR_FN}"; }
fi

echo ""
echo "=== Download summary ==="
ls -lh "${DEST}/"*.tsv.bgz 2>&1 | head -20
echo ""
echo "End: $(date)"
