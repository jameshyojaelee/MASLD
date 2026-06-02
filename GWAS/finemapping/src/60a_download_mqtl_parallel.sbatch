#!/bin/bash
#SBATCH --job-name=mqtldl
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --array=0-15
#SBATCH --output=GWAS/finemapping/logs/60a_dl_par_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/60a_dl_par_%A_%a.err
# ---------------------------------------------------------------------------
# 60a_download_mqtl_parallel.sbatch
# Parallel (8-way strided) replacement for the serial 60a downloader.
# Each array task downloads its 1/8 stride of the combined Chen+Ottensmann list.
# Skip-if-exists with gzip-integrity check (idempotent — keeps the ~29 files the
# serial job already fetched, re-fetches any partial/corrupt). Atomic via curl -o
# + gzip -t verification before keeping. 8 streams -> ~5-8x the serial rate
# (bounded by EBI per-IP limits; 8 is the polite ceiling).
# ---------------------------------------------------------------------------
set -uo pipefail

BASE_DIR="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$BASE_DIR" || exit 1

CHEN_DIR="$BASE_DIR/data/external/chen2023_mqtl/sumstats"
LIPID_DIR="$BASE_DIR/data/external/lipidqtl/sumstats"
CHEN_LIST="$BASE_DIR/data/external/chen2023_mqtl/chen2023_biomarker_download_list.tsv"
LIPID_LIST="$BASE_DIR/data/external/lipidqtl/ottensmann2023_accession_trait_map.tsv"
mkdir -p "$CHEN_DIR" "$LIPID_DIR" "$BASE_DIR/GWAS/finemapping/logs"

FTP_ROOT="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
TASK="${SLURM_ARRAY_TASK_ID:-0}"
NTASKS="${SLURM_ARRAY_TASK_COUNT:-8}"

acc_block() {
  local acc="$1"; local num="${acc#GCST}"; num=$((10#$num))
  local lo=$(( (num - 1) / 1000 * 1000 + 1 )); local hi=$(( lo + 999 ))
  printf "GCST%08d-GCST%08d" "$lo" "$hi"
}
download_one() {
  local acc="$1" dest="$2"
  local block; block=$(acc_block "$acc")
  local base_url="$FTP_ROOT/$block/$acc"
  local candidates=("${acc}_buildGRCh38.tsv.gz" "${acc}.tsv.gz" "${acc}_buildGRCh37.tsv.gz")
  local outf="$dest/${acc}.tsv.gz"
  if [[ -s "$outf" ]] && gzip -t "$outf" 2>/dev/null; then echo "  [skip] $acc"; return 0; fi
  [[ -e "$outf" ]] && rm -f "$outf"
  for fn in "${candidates[@]}"; do
    local url="$base_url/$fn"
    if curl -sf -I "$url" >/dev/null 2>&1; then
      echo "  [get ] task$TASK $acc <- $fn"
      if curl -sf --retry 4 --retry-delay 5 -o "$outf" "$url" && gzip -t "$outf" 2>/dev/null; then
        return 0
      fi
      [[ -e "$outf" ]] && rm -f "$outf"
    fi
  done
  echo "  [FAIL] $acc (no candidate downloaded)"; return 1
}

# Unified work list: "<accession>\t<dest_dir>"  (Chen then Ottensmann)
WORK=$(mktemp)
{ tail -n +2 "$CHEN_LIST"  | awk -v d="$CHEN_DIR"  'NF{print $1"\t"d}';
  tail -n +2 "$LIPID_LIST" | awk -v d="$LIPID_DIR" 'NF{print $1"\t"d}'; } > "$WORK"
total=$(wc -l < "$WORK")
echo "task $TASK/$NTASKS : $total accessions total; processing stride"

ok=0; fail=0; i=0
while IFS=$'\t' read -r acc dest; do
  if (( i % NTASKS == TASK )); then
    if download_one "$acc" "$dest"; then ok=$((ok+1)); else fail=$((fail+1)); fi
  fi
  i=$((i+1))
done < "$WORK"
rm -f "$WORK"
echo "task $TASK done: ok=$ok fail=$fail"
