#!/bin/bash
# 77_build_eve_merged.sh -- extract per-protein EVE VCFs from the evemodel.org
# bulk zip and merge into ONE bgzip+tabix GRCh38 VCF for coordinate lookup.
# EVE per-protein VCFs are UniProt-entry-named, so we join by coord not gene.
# Additive; feeds 77_coding_upgrade.py (--eve-merged).
set -euo pipefail
ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
EVE="$ROOT/data/external/eve"
ZIP="$EVE/eve_bulk.zip"
EXTRACT="$EVE/extracted"
MERGED="$EVE/eve_merged.vcf"
export PATH="/nfs/sw/easybuild/software/htslib/1.23.1/bin:$PATH"

echo "[1/5] verifying zip"
unzip -t "$ZIP" >/dev/null

echo "[2/5] locating VCF folder inside archive"
VCF_DIR=$(unzip -Z1 "$ZIP" | grep -iE 'vcf.*\.vcf$' | head -1 | xargs -r dirname)
[ -z "$VCF_DIR" ] && { echo "no VCF folder found in zip"; exit 1; }
echo "    -> $VCF_DIR"

echo "[3/5] extracting only $VCF_DIR/*.vcf"
rm -rf "$EXTRACT"; mkdir -p "$EXTRACT"
unzip -o -j "$ZIP" "$VCF_DIR/*.vcf" -d "$EXTRACT" >/dev/null
N=$(ls "$EXTRACT"/*.vcf 2>/dev/null | wc -l)
echo "    -> $N per-protein VCFs"
[ "$N" -eq 0 ] && { echo "no VCFs extracted"; exit 1; }

echo "[4/5] merging (header from first file; chrom normalised to no-chr)"
FIRST=$(ls "$EXTRACT"/*.vcf | head -1)
grep '^#' "$FIRST" > "$MERGED"
# body: strip 'chr', sort by chrom/pos; chrom may be 1..22,X,Y
cat "$EXTRACT"/*.vcf | grep -v '^#' \
  | sed 's/^chr//' \
  | sort -t$'\t' -k1,1 -k2,2n >> "$MERGED"
BODY=$(grep -vc '^#' "$MERGED")
echo "    -> $BODY total EVE variant rows"

echo "[5/5] bgzip + tabix"
bgzip -f "$MERGED"
tabix -f -p vcf "$MERGED.gz"
echo "DONE: $MERGED.gz (+ .tbi)"
