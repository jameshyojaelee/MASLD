#!/usr/bin/env bash
#SBATCH --job-name=wget
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/refetch_hmsma_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/refetch_hmsma_%j.out
##############################################################################
# Re-fetch the 5 HMSMA FASTQs whose gzip stream is TRUNCATED.
# Verified 2026-08-09 by `gzip -t` + read counting (_integrity_report.tsv):
# each has a corrupt member and a read count far below its mate, which is why
# STARsolo aborted with FASTQ record desynchronisation.
#
#   HRA_19  HRR1782845_r2   27.3M reads vs mate 226.5M
#   HRA_24  HRR1782853_r2   68.6M       vs mate 291.9M
#   HRA_27  HRR2173370_f1  153.0M       vs mate 1232.1M
#   HRA_28  HRR2173371_r2  156.7M       vs mate 445.0M
#   HRA_50  HRR2173358_r2  146.5M       vs mate 285.6M
#
# Corrupt originals are RENAMED to *.truncated_<ts>, never deleted.
# Re-verifies gzip integrity after download and reports per file.
##############################################################################
set -o pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUTDIR="data/HRA007511/fastq_visium"
FTP_BASE="ftp://download.big.ac.cn/gsa-human/HRA007511"
TS=$(date +%Y%m%d_%H%M%S)

# run_id:file_suffix
TARGETS=(
  "HRR1782845:r2"
  "HRR1782853:r2"
  "HRR2173370:f1"
  "HRR2173371:r2"
  "HRR2173358:r2"
)

FAIL=0
for T in "${TARGETS[@]}"; do
  RUN="${T%%:*}"; SUF="${T##*:}"
  F="${OUTDIR}/${RUN}/${RUN}_${SUF}.fastq.gz"
  URL="${FTP_BASE}/${RUN}/${RUN}_${SUF}.fastq.gz"

  if [ -f "$F" ]; then
    mv "$F" "${F}.truncated_${TS}"
    echo "ARCHIVED  ${F} -> $(basename "${F}").truncated_${TS}"
  fi

  echo "FETCHING  $URL"
  # NB: no --continue; it conflicts with -O and the target was just archived.
  wget -q --tries=5 --retry-connrefused --timeout=60 \
       -O "$F" "$URL"
  RC=$?
  if [ $RC -ne 0 ]; then
    echo "DOWNLOAD_FAILED ${RUN}_${SUF} rc=$RC"; FAIL=$((FAIL+1)); continue
  fi

  if gzip -t "$F" 2>/dev/null; then
    N=$(zcat "$F" 2>/dev/null | wc -l); R=$((N/4))
    echo "OK        ${RUN}_${SUF}  reads=${R}  size=$(du -h "$F" | cut -f1)"
  else
    echo "STILL_CORRUPT ${RUN}_${SUF} — source file may itself be bad"; FAIL=$((FAIL+1))
  fi
done

echo "REFETCH_DONE failures=${FAIL}"
[ "$FAIL" -eq 0 ] && echo "READY_FOR_STARSOLO_RERUN indices=8,13,15,16,23"
