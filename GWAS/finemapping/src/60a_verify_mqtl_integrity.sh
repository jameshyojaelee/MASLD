#!/bin/bash
#SBATCH --job-name=mqtlqc
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=GWAS/finemapping/logs/60a_verify_%j.out
#SBATCH --error=GWAS/finemapping/logs/60a_verify_%j.err
# ---------------------------------------------------------------------------
# 60a_verify_mqtl_integrity.sh
# Integrity gate between the download array and the COLOC arrays.
# gzip -t every expected file (catches truncation from the cancel/resume of the
# serial downloader); re-download any corrupt/missing; sanity-check >1 data row.
# Exits NON-ZERO if any file is still unrepairable -> the afterok COLOC arrays
# will NOT fire on bad data.
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
  done
  return 1
}

WORK=$(mktemp)
{ tail -n +2 "$CHEN_LIST"  | awk -v d="$CHEN_DIR"  'NF{print $1"\t"d}';
  tail -n +2 "$LIPID_LIST" | awk -v d="$LIPID_DIR" 'NF{print $1"\t"d}'; } > "$WORK"

ok=0; corrupt=0; missing=0; repaired=0; unrepairable=0; emptyrows=0
while IFS=$'\t' read -r acc dest; do
  outf="$dest/${acc}.tsv.gz"
  if [[ ! -s "$outf" ]]; then
    missing=$((missing+1))
    if download_one "$acc" "$dest"; then repaired=$((repaired+1)); else unrepairable=$((unrepairable+1)); echo "  [UNREPAIRABLE-missing] $acc"; continue; fi
  elif ! gzip -t "$outf" 2>/dev/null; then
    corrupt=$((corrupt+1)); echo "  [corrupt] $acc -> re-download"; rm -f "$outf"
    if download_one "$acc" "$dest"; then repaired=$((repaired+1)); else unrepairable=$((unrepairable+1)); echo "  [UNREPAIRABLE-corrupt] $acc"; continue; fi
  else
    ok=$((ok+1))
  fi
  # sanity: >1 line (header + >=1 data row)
  nl=$(zcat "$outf" 2>/dev/null | head -3 | wc -l)
  if [[ "$nl" -lt 2 ]]; then emptyrows=$((emptyrows+1)); echo "  [EMPTY] $acc (<2 lines)"; fi
done < "$WORK"
rm -f "$WORK"

echo "============================================================"
echo "INTEGRITY: ok=$ok corrupt=$corrupt missing=$missing repaired=$repaired unrepairable=$unrepairable empty=$emptyrows"
echo "Chen files:  $(ls "$CHEN_DIR"/*.tsv.gz 2>/dev/null | wc -l)"
echo "Lipid files: $(ls "$LIPID_DIR"/*.tsv.gz 2>/dev/null | wc -l)"
echo "============================================================"
# Fail the gate if anything is still bad -> COLOC (afterok) will not run on bad data
if [[ "$unrepairable" -gt 0 || "$emptyrows" -gt 0 ]]; then
  echo "GATE FAILED: $unrepairable unrepairable + $emptyrows empty"; exit 1
fi
echo "GATE PASSED: all files gzip-valid with data."
