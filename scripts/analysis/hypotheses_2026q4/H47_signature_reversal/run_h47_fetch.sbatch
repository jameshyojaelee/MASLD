#!/bin/bash
#SBATCH --job-name=H47fetch
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
# H47 step A: download LINCS level 5 GCTX (GSE92742, GSE70138), check GEO SHA512,
# decompress, map landmarks to refit Ensembl IDs, extract landmark x liver-line
# compound signatures. No scoring.
# Usage: sbatch --output=<H47>/logs/slurm-%j.log run_h47_fetch.sbatch <H47_dir> <extract_dir>
set -euo pipefail
H47="${1:?H47 dir required}"
OUT="${2:?extract dir required}"
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H47_signature_reversal
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
DEG=$ROOT/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv
META=$ROOT/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz
DL="$H47/downloads_untrusted"
DEC="$H47/downloads_untrusted_decompressed"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS_USER= PYTHONHASHSEED=0
module load pigz/2.7-GCCcore-12.2.0 >/dev/null 2>&1 || true
test -s "$H47/PRESPEC_H47.md"
test ! -e "$OUT"
mkdir -p "$DEC"

B1=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl
B2=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl
F1=GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz
F2=GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz

cd "$DL"
zcat GSE92742_SHA512SUMS.txt.gz GSE70138_SHA512SUMS.txt.gz | grep -E "  ($F1|$F2)\$" > level5_SHA512SUMS.txt
test "$(wc -l < level5_SHA512SUMS.txt)" -eq 2
for pair in "$B1 $F1" "$B2 $F2"; do
  set -- $pair
  if [ ! -s "$2" ] || ! grep "  $2\$" level5_SHA512SUMS.txt | sha512sum -c --quiet - ; then
    echo "[$(date -u +%FT%TZ)] downloading $2"
    curl -sS --fail --retry 5 --retry-delay 30 -C - -o "$2" "$1/$2"
  fi
done
echo "[$(date -u +%FT%TZ)] SHA512 check (stop rule S1)"
sha512sum -c level5_SHA512SUMS.txt

for f in $F1 $F2; do
  g="$DEC/${f%.gz}"
  if [ ! -s "$g" ]; then
    echo "[$(date -u +%FT%TZ)] decompressing $f"
    pigz -dc "$f" > "$g.partial" 2>/dev/null || gzip -dc "$f" > "$g.partial"
    mv "$g.partial" "$g"
  fi
done
cd "$DEC"
sha256sum ./*.gctx > "$DEC/decompressed_sha256.txt"

echo "[$(date -u +%FT%TZ)] extracting"
"$ENV/bin/python3" -I "$SCRIPTS/h47_extract.py" "$DL" "$DEC" "$OUT"
echo "[$(date -u +%FT%TZ)] mapping landmarks"
"$ENV/bin/Rscript" --vanilla "$SCRIPTS/h47_map_landmarks.R" "$DL/GSE92742_Broad_LINCS_gene_info.txt.gz" "$DEG" "$META" "$OUT"
(cd "$SCRIPTS" && sha256sum h47_extract.py h47_map_landmarks.R run_h47_fetch.sbatch) > "$OUT/script_sha256.txt"
sha256sum "$DL"/*.gz "$DL"/*.txt "$DEC/decompressed_sha256.txt" > "$OUT/input_sha256.txt"
cp "$DEC/decompressed_sha256.txt" "$OUT/"
"$ENV/bin/python3" -I -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
echo "[$(date -u +%FT%TZ)] done"
