#!/bin/bash
#SBATCH --job-name=finngen_dl
#SBATCH --partition=cpu,io
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/finngen_dl_%A_%a.out
#SBATCH --error=logs/finngen_dl_%A_%a.err

set -euo pipefail

TRAITS=(NAFLD CHIRHEP_NAS C3_HEPATOCELLU_CARC_EXALLC E4_OBESITY)
TRAIT=${TRAITS[$SLURM_ARRAY_TASK_ID]}

OUTDIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/FinnGen
mkdir -p "$OUTDIR"

URL="https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_${TRAIT}.gz"
OUTFILE="${OUTDIR}/finngen_R12_${TRAIT}.gz"

echo "=== FinnGen Download: ${TRAIT} ==="
echo "Job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Partition: ${SLURM_JOB_PARTITION}"
echo "Start: $(date)"
echo "URL: ${URL}"

wget -c -q --show-progress -O "$OUTFILE" "$URL"

# Verify file is non-empty and gzip-valid
if [ ! -s "$OUTFILE" ]; then
  echo "ERROR: Downloaded file is empty"
  exit 1
fi

gzip -t "$OUTFILE" 2>/dev/null
if [ $? -ne 0 ]; then
  echo "ERROR: File failed gzip integrity check"
  exit 1
fi

echo "OK: $(ls -lh "$OUTFILE")"
echo "End: $(date)"
