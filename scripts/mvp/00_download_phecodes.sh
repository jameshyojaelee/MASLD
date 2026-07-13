#!/bin/bash
#SBATCH --job-name=wget
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/logs/download_phecodes_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/logs/download_phecodes_%j.err

# 00_download_phecodes.sh
# Download the MVP R4 PheCodes DigestiveSystem_batch2 tar (200.7 GB, dbGaP phs002453,
# open access) and extract the liver phecode (571.x) sumstats + ALL metadata.
# Aspera auth is dead on this cluster -> wget -c (resumable). md5-verified.

set -euo pipefail

BASE_URL="https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/phs002453.v1.p1/analyses/GIA"
TARNAME="phs002453.MVP_R4.1000G_AGR.GIA.PheCodes_DigestiveSystem_batch2.analysis-PI.MULTI.tar"
EXPECTED_MD5="6a5aafa8ab2a2447aee78f6c3ef02655"

SCRATCH="/scratch/claude-91825/-gpfs-commons-groups-sanjana-lab-Cas13-MASLD-library-design/aeeb258d-a882-4693-b330-1abbd61c1b7d/mvp_dl"
DEST="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/raw/phecodes"
mkdir -p "$SCRATCH" "$DEST"
cd "$SCRATCH"

echo "[$(date)] Downloading $TARNAME (~200 GB) via wget -c ..."
time wget -c -q "$BASE_URL/$TARNAME" -O "$TARNAME"
echo "[$(date)] Download finished. Size: $(ls -lh $TARNAME | awk '{print $5}')"

echo "[$(date)] Verifying md5 ..."
GOT_MD5=$(md5sum "$TARNAME" | awk '{print $1}')
if [ "$GOT_MD5" != "$EXPECTED_MD5" ]; then
    echo "ERROR: md5 mismatch! got=$GOT_MD5 expected=$EXPECTED_MD5"
    exit 1
fi
echo "  md5 OK ($GOT_MD5)"

echo "[$(date)] Extracting Phe_571* (sumstats + metadata) + all metadata ..."
# Extract liver phecode 571.x gz + their metadata (and any 573 metadata for completeness)
tar -xvf "$TARNAME" -C "$DEST" --strip-components=1 --wildcards \
    '*Phe_571*' '*Phe_573*.metadata.txt' 2>&1 | tail -40

echo "[$(date)] Integrity-checking extracted .gz files ..."
NFAIL=0
for f in "$DEST"/*.txt.gz; do
    if ! gzip -t "$f" 2>/dev/null; then echo "  CORRUPT: $f"; NFAIL=$((NFAIL+1)); fi
done
echo "  gz integrity: $NFAIL failures"
echo "  extracted files: $(ls "$DEST" | wc -l)"
[ "$NFAIL" -eq 0 ] || { echo "ERROR: corrupt gz extracted"; exit 1; }

echo "[$(date)] Deleting 200 GB tar to free /scratch ..."
rm -f "$TARNAME"

echo "[$(date)] DONE. Liver phecode files staged to $DEST"
ls -lh "$DEST" | head -50
