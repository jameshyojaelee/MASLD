#!/bin/bash
#SBATCH --job-name=mqtlqc
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=GWAS/finemapping/logs/60a_verify_par_%j.out
#SBATCH --error=GWAS/finemapping/logs/60a_verify_par_%j.err
# ---------------------------------------------------------------------------
# 60a_verify_mqtl_integrity_parallel.sh
# Parallel (16-way) integrity gate — same guarantees as the serial version but
# gzip -t runs concurrently across files. Files already downloaded; this only
# tests + repairs. Exits NON-ZERO if anything is unrepairable so the afterok
# COLOC arrays will not run on bad data.
# ---------------------------------------------------------------------------
set -uo pipefail
BASE_DIR="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$BASE_DIR" || exit 1
CHEN_DIR="$BASE_DIR/data/external/chen2023_mqtl/sumstats"
LIPID_DIR="$BASE_DIR/data/external/lipidqtl/sumstats"
CHEN_LIST="$BASE_DIR/data/external/chen2023_mqtl/chen2023_biomarker_download_list.tsv"
LIPID_LIST="$BASE_DIR/data/external/lipidqtl/ottensmann2023_accession_trait_map.tsv"
FTP_ROOT="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

acc_block() { local acc="$1"; local num="${acc#GCST}"; num=$((10#$num));
  local lo=$(( (num-1)/1000*1000+1 )); printf "GCST%08d-GCST%08d" "$lo" $((lo+999)); }
download_one() {
  local acc="$1" dest="$2"; local block; block=$(acc_block "$acc")
  local base_url="$FTP_ROOT/$block/$acc" outf="$dest/${acc}.tsv.gz"
  for fn in "${acc}_buildGRCh38.tsv.gz" "${acc}.tsv.gz" "${acc}_buildGRCh37.tsv.gz"; do
    if curl -sf -I "$base_url/$fn" >/dev/null 2>&1; then
      curl -sf --retry 5 --retry-delay 5 -o "$outf" "$base_url/$fn" && gzip -t "$outf" 2>/dev/null && return 0
      [[ -e "$outf" ]] && rm -f "$outf"
    fi
  done; return 1
}

WORK=$(mktemp)
{ tail -n +2 "$CHEN_LIST"  | awk -v d="$CHEN_DIR"  'NF{print $1"\t"d}';
  tail -n +2 "$LIPID_LIST" | awk -v d="$LIPID_DIR" 'NF{print $1"\t"d}'; } > "$WORK"
ntot=$(wc -l < "$WORK")

# ---- Phase 1: parallel gzip -t (the bottleneck), 16-way ----
PROB=$(mktemp)
test_one() {
  local acc="$1" dest="$2" f="$2/$1.tsv.gz"
  if [[ ! -s "$f" ]]; then echo "MISSING $acc $dest"; return; fi
  if ! gzip -t "$f" 2>/dev/null; then echo "CORRUPT $acc $dest"; return; fi
  if [[ $(zcat "$f" 2>/dev/null | head -2 | wc -l) -lt 2 ]]; then echo "EMPTY $acc $dest"; fi
}
export -f test_one
awk -F'\t' '{print $1" "$2}' "$WORK" | xargs -P 16 -n 2 bash -c 'test_one "$@"' _ > "$PROB" 2>/dev/null

nprob=$(grep -c . "$PROB" 2>/dev/null || true); nprob=${nprob:-0}
ok=$(( ntot - nprob ))
echo "Phase 1 (parallel gzip -t): $ntot files, $ok ok, $nprob flagged"

# ---- Phase 2: repair flagged (expected ~0; sequential, few) ----
repaired=0; unrepairable=0; empty=0
while read -r status acc dest; do
  [[ -z "${status:-}" ]] && continue
  case "$status" in
    EMPTY) empty=$((empty+1)); echo "  [EMPTY] $acc";;
    MISSING|CORRUPT)
      echo "  [$status] $acc -> re-download"; rm -f "$dest/$acc.tsv.gz"
      if download_one "$acc" "$dest"; then repaired=$((repaired+1)); else unrepairable=$((unrepairable+1)); echo "  [UNREPAIRABLE] $acc"; fi;;
  esac
done < "$PROB"
rm -f "$WORK" "$PROB"

echo "============================================================"
echo "INTEGRITY: total=$ntot ok=$ok flagged=$nprob repaired=$repaired unrepairable=$unrepairable empty=$empty"
echo "Chen files:  $(ls "$CHEN_DIR"/*.tsv.gz 2>/dev/null | wc -l)"
echo "Lipid files: $(ls "$LIPID_DIR"/*.tsv.gz 2>/dev/null | wc -l)"
echo "============================================================"
if [[ "$unrepairable" -gt 0 || "$empty" -gt 0 ]]; then
  echo "GATE FAILED: $unrepairable unrepairable + $empty empty"; exit 1
fi
echo "GATE PASSED: all files gzip-valid with data."
