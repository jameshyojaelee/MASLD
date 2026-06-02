#!/bin/bash
#SBATCH --job-name=coloc
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=GWAS/finemapping/logs/60a_download_mqtl_%j.out
#SBATCH --error=GWAS/finemapping/logs/60a_download_mqtl_%j.err
# ---------------------------------------------------------------------------
# 60a_download_mqtl_sumstats.sh
# Download metabolite-QTL (Chen 2023, PMID 36635386) + lipid-species-QTL
# (Ottensmann 2023, PMID 37907536) GWAS summary stats from EBI GWAS Catalog FTP.
#
# Both are GRCh38 (hg38). Liftover to hg19 happens inside 60b/60c per-locus
# (the MASLD GWAS sumstats are hg19).
#
# Chen: restrict to MASLD-biomarker metabolites only (BCAAs / ceramides /
#       bile acids / acylcarnitines / PC-PE) to bound the download.
#       Download list: data/external/chen2023_mqtl/chen2023_biomarker_download_list.tsv
# Ottensmann: all 179 lipid species (the entire layer IS lipid species).
#
# This is an I/O-heavy job -> io partition, NOT the login node.
# ---------------------------------------------------------------------------
set -uo pipefail

BASE_DIR="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$BASE_DIR" || exit 1

CHEN_DIR="$BASE_DIR/data/external/chen2023_mqtl/sumstats"
LIPID_DIR="$BASE_DIR/data/external/lipidqtl/sumstats"
mkdir -p "$CHEN_DIR" "$LIPID_DIR"

FTP_ROOT="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

# GWAS Catalog FTP groups accessions in blocks of 1000:
#   GCST90199621 -> GCST90199001-GCST90200000
acc_block() {
  local acc="$1"
  local num="${acc#GCST}"
  num=$((10#$num))
  local lo=$(( (num - 1) / 1000 * 1000 + 1 ))
  local hi=$(( lo + 999 ))
  printf "GCST%08d-GCST%08d" "$lo" "$hi"
}

download_one() {
  # $1 accession, $2 dest dir, $3 build-suffix-pattern (file name guesses tried in order)
  local acc="$1" dest="$2"
  local block; block=$(acc_block "$acc")
  local base_url="$FTP_ROOT/$block/$acc"
  # Candidate file names (Chen uses *_buildGRCh38.tsv.gz; Ottensmann uses ACC.tsv.gz)
  local candidates=("${acc}_buildGRCh38.tsv.gz" "${acc}.tsv.gz" "${acc}_buildGRCh37.tsv.gz")
  local outf="$dest/${acc}.tsv.gz"
  if [[ -s "$outf" ]]; then
    # verify gzip integrity; if good, skip
    if gzip -t "$outf" 2>/dev/null; then echo "  [skip] $acc (already present)"; return 0; fi
    echo "  [redo] $acc (corrupt) "; rm -f "$outf"
  fi
  for fn in "${candidates[@]}"; do
    local url="$base_url/$fn"
    # HEAD check
    if curl -sf -I "$url" >/dev/null 2>&1; then
      echo "  [get ] $acc <- $fn"
      if curl -sf --retry 4 --retry-delay 5 -o "$outf" "$url"; then
        if gzip -t "$outf" 2>/dev/null; then return 0; fi
        echo "  [warn] $acc downloaded but failed gzip -t; removing"; rm -f "$outf"
      fi
    fi
  done
  echo "  [FAIL] $acc : no candidate file downloaded"
  return 1
}

echo "============================================================"
echo "60a metabolite/lipid QTL download"
echo "start: $(date)"
echo "============================================================"

# --- Chen 2023 biomarker metabolites ---
CHEN_LIST="$BASE_DIR/data/external/chen2023_mqtl/chen2023_biomarker_download_list.tsv"
if [[ ! -s "$CHEN_LIST" ]]; then
  echo "ERROR: Chen download list missing: $CHEN_LIST"; exit 1
fi
echo ""
echo "--- Chen 2023 metabolites ($(($(wc -l < "$CHEN_LIST")-1)) biomarker traits) ---"
chen_ok=0; chen_fail=0
while IFS=$'\t' read -r acc trait nnum; do
  [[ "$acc" == "accession" ]] && continue
  if download_one "$acc" "$CHEN_DIR"; then chen_ok=$((chen_ok+1)); else chen_fail=$((chen_fail+1)); fi
done < "$CHEN_LIST"
echo "Chen done: ok=$chen_ok fail=$chen_fail"

# --- Ottensmann 2023 lipid species ---
LIPID_LIST="$BASE_DIR/data/external/lipidqtl/ottensmann2023_accession_trait_map.tsv"
if [[ ! -s "$LIPID_LIST" ]]; then
  echo "ERROR: Ottensmann map missing: $LIPID_LIST"; exit 1
fi
echo ""
echo "--- Ottensmann 2023 lipid species ($(($(wc -l < "$LIPID_LIST")-1)) traits) ---"
lip_ok=0; lip_fail=0
while IFS=$'\t' read -r acc trait sample; do
  [[ "$acc" == "accession" ]] && continue
  if download_one "$acc" "$LIPID_DIR"; then lip_ok=$((lip_ok+1)); else lip_fail=$((lip_fail+1)); fi
done < "$LIPID_LIST"
echo "Ottensmann done: ok=$lip_ok fail=$lip_fail"

echo ""
echo "============================================================"
echo "Chen files:        $(ls "$CHEN_DIR"/*.tsv.gz 2>/dev/null | wc -l)"
echo "Ottensmann files:  $(ls "$LIPID_DIR"/*.tsv.gz 2>/dev/null | wc -l)"
echo "Chen disk:         $(du -sh "$CHEN_DIR" 2>/dev/null | cut -f1)"
echo "Ottensmann disk:   $(du -sh "$LIPID_DIR" 2>/dev/null | cut -f1)"
echo "end: $(date)"
echo "============================================================"
# Touch a sentinel so the COLOC array knows the download finished
touch "$BASE_DIR/data/external/chen2023_mqtl/.download_complete"
