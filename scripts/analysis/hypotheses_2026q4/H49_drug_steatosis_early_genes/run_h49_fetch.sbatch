#!/bin/bash
#SBATCH --job-name=H49fetch
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=90:00:00
# H49 step A: download the 158 Open TG-GATEs human in vitro CEL archives,
# 8 concurrent connections, study drugs and controls first.
# Usage: sbatch --output=<H49>/logs/fetch-%j.log run_h49_fetch.sbatch <H49_dir>
set -euo pipefail
H49="${1:?H49 dir required}"
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H49_drug_steatosis_early_genes
DL="$H49/downloads_untrusted_tggates_cel"
test -s "$H49/PRESPEC_H49.md"
mkdir -p "$DL"
URLS="$H49/logs/tggates_human_urls.tsv"
[ -s "$URLS" ] || "$ENV/bin/python3" -I "$SCRIPTS/h49_url_list.py" \
  "$H49/downloads_untrusted_tggates_meta_unzipped/open_tggates_main.csv" "$URLS"
PRIORITY='^(amiodarone|tamoxifen|perhexiline|tetracycline|nefazodone|carbamazepine|enalapril|isoniazid|ketoconazole|valproic_acid|omeprazole|caffeine|chloramphenicol|chlorpheniramine|hydroxyzine|theophylline|rifampicin|acetaminophen|ibuprofen|naproxen|aspirin|dexamethasone)\.'
ORDERED="$H49/logs/tggates_human_urls_ordered.tsv"
{ grep -E "$PRIORITY" "$URLS" || true; grep -vE "$PRIORITY" "$URLS" || true; } > "$ORDERED"
test "$(wc -l < "$ORDERED")" = 158

fetch_one() {
  local name="$1" url="$2" dest="$DL/$1"
  if [ -s "$dest" ] && unzip -tq "$dest" >/dev/null 2>&1; then echo "have $name"; return 0; fi
  for attempt in 1 2 3 4 5; do
    if curl -sS -L --fail -m 7200 -o "$dest.part" "$url" && unzip -tq "$dest.part" >/dev/null 2>&1; then
      mv "$dest.part" "$dest"; echo "ok $name"; return 0
    fi
    echo "retry $attempt $name" >&2; sleep 30
  done
  echo "FAILED $name" >&2; return 1
}
export -f fetch_one
export DL
set +e
awk -F'\t' '{print $1"\t"$2}' "$ORDERED" | tr '\t' '\n' | xargs -d '\n' -n 2 -P 8 bash -c 'fetch_one "$0" "$1"'
rc=$?
set -e
n_ok=$(ls "$DL"/*.zip 2>/dev/null | wc -l)
echo "archives present: $n_ok of 158 (xargs rc=$rc)"
(cd "$DL" && sha256sum *.zip > "$H49/logs/tggates_cel_zip_sha256.txt")
du -sh "$DL"
test "$n_ok" = 158
