#!/bin/bash
#SBATCH --job-name=wget
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=6
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/logs/download_labs_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/logs/download_labs_%j.err

# 00_download_labs.sh
# Download MVP R4 liver-biomarker GWAS (Mean_INT) individually from the expanded
# dbGaP Labs dirs: ALT/AST/Albumin (batch1) + Platelet (batch3), x EUR/AFR/AMR/EAS/META.
# 20 sumstat .gz + 20 metadata.txt. (ALT+AST+Platelet = FIB-4 fibrosis components.)

set -euo pipefail

BASE="https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/phs002453.v1.p1/analyses/GIA"
DEST="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/raw/labs"
mkdir -p "$DEST"

# analyte -> Labs batch
declare -A BATCH=( [ALT]=1 [AST]=1 [Albumin]=1 [Platelet]=3 )

# Build the URL list
URLS=$(mktemp)
for analyte in ALT AST Albumin Platelet; do
  b=${BATCH[$analyte]}
  LABS="$BASE/MVP_R4.1000G_AGR.GIA.Labs_batch${b}.tar_expanded/Labs"
  for anc in EUR AFR AMR EAS META; do
    for ext in txt.gz metadata.txt; do
      echo "$LABS/MVP_R4.1000G_AGR.${analyte}_Mean_INT.${anc}.GIA.dbGaP.${ext}"
    done
  done
done > "$URLS"

echo "[$(date)] Downloading $(wc -l < "$URLS") MVP liver-biomarker files (parallel x5) ..."
cd "$DEST"
# -P5 parallel; -c resume; -q quiet
cat "$URLS" | xargs -P 5 -I{} wget -c -q {}

echo "[$(date)] Integrity-checking .gz ..."
NFAIL=0
for f in "$DEST"/*.txt.gz; do
  if ! gzip -t "$f" 2>/dev/null; then echo "  CORRUPT: $f"; NFAIL=$((NFAIL+1)); fi
done
echo "  files: $(ls "$DEST" | wc -l) | gz: $(ls "$DEST"/*.txt.gz 2>/dev/null | wc -l) | corrupt: $NFAIL"
rm -f "$URLS"
[ "$NFAIL" -eq 0 ] || { echo "ERROR: corrupt gz"; exit 1; }
echo "[$(date)] DONE. Labs staged to $DEST"
ls -lh "$DEST"/*.txt.gz | awk '{print $5, $9}'
