#!/bin/bash
#SBATCH --job-name=V3dl
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
# Download the 17 MASLD-spectrum human GSE256398 CellBender h5 files and the
# GEO family SOFT metadata into an empty untrusted directory.
# Usage: sbatch download_gse256398.sbatch <DEST_DIR>
set -euo pipefail
DEST="${1:?dest dir required}"
if [ -n "$(ls -A "$DEST")" ]; then echo "destination not empty" >&2; exit 1; fi
cd "$DEST"
FILES="GSM8097080_S12 GSM8097081_S13 GSM8097082_S14 GSM8097083_S15 GSM8097084_S16
GSM8097085_S17 GSM8097086_S18 GSM8097087_S19 GSM8097088_S20 GSM8097089_S21
GSM8097090_S23 GSM8097091_S27 GSM8097092_S28 GSM8097093_S33 GSM8097094_S35
GSM8097095_S36 GSM8097096_S38"
for f in $FILES; do
  gsm=${f%%_*}
  url="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8097nnn/${gsm}/suppl/${f}_CB_raw_feature_bc_matrix_filtered.h5"
  curl -sS --fail --retry 5 --retry-delay 20 -o "${f}_CB_raw_feature_bc_matrix_filtered.h5" "$url"
done
curl -sS --fail --retry 5 -o GSE256398_family.soft.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256398/soft/GSE256398_family.soft.gz"
curl -sS --fail --retry 5 -o filelist.txt \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256398/suppl/filelist.txt"
sha256sum *.h5 GSE256398_family.soft.gz > SHA256SUMS
ls -l
