#!/usr/bin/env bash
#SBATCH --job-name=wget
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/diag_hrr2173371_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/diag_hrr2173371_%j.out
##############################################################################
# HRR2173371_r2 (HRA_28 R2) re-fetched cleanly on read count — 444,952,267
# reads, exactly matching its R1 — but still fails `gzip -t`. Data appears
# COMPLETE; the gzip stream trailer is bad. Distinguish transfer artefact from
# a genuinely malformed source file by capturing the exact error, re-fetching
# once, and re-testing.
##############################################################################
set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
D=data/HRA007511/fastq_visium/HRR2173371
F=$D/HRR2173371_r2.fastq.gz
URL="ftp://download.big.ac.cn/gsa-human/HRA007511/HRR2173371/HRR2173371_r2.fastq.gz"
TS=$(date +%Y%m%d_%H%M%S)

echo "=== attempt 1: exact gzip error on the current file ==="
gzip -t "$F" 2>&1 | head -5
echo "  size: $(stat -c %s "$F") bytes"

echo "=== re-fetch once ==="
mv "$F" "${F}.badtrailer_${TS}"
wget -q --tries=5 --retry-connrefused --timeout=60 -O "$F" "$URL"
echo "  wget rc=$?  size: $(stat -c %s "$F") bytes"

echo "=== attempt 2: gzip test after re-fetch ==="
if gzip -t "$F" 2>gz_err.txt; then
  echo "GZIP_OK_AFTER_REFETCH"
  VERDICT=transfer_artefact
else
  echo "GZIP_STILL_BAD: $(head -2 gz_err.txt)"
  VERDICT=source_file_malformed
fi
rm -f gz_err.txt

N=$(zcat "$F" 2>/dev/null | wc -l); R=$((N/4))
echo "  reads after re-fetch: $R (R1 mate has 444952267)"
[ "$R" -eq 444952267 ] && echo "  READ_COUNT_MATCHES_MATE" || echo "  READ_COUNT_MISMATCH"

echo "VERDICT=${VERDICT}"
if [ "$VERDICT" = "source_file_malformed" ] && [ "$R" -eq 444952267 ]; then
  echo "GUIDANCE: data is complete and pair-matched; the defect is confined to the"
  echo "gzip trailer. STAR reads via zcat and will consume all records. Proceeding"
  echo "is defensible, but must be recorded as a known source-file defect."
fi
