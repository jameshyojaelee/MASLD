#!/bin/bash
# Parallel array download of the full Sveinbjornsson-2022 NAFLD bundle (decode.com data-use form).
# One file per array task; md5-verified; skips files already present with a matching checksum.
#SBATCH --job-name=decode_dl
#SBATCH --partition=io
#SBATCH --array=1-12%8
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --output=logs/decode_dl_%A_%a.out
set -uo pipefail
DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/Sveinbjornsson2022
cd "$DIR"
TOKEN=06828a13-b8b6-46b7-bfe1-62d1b4eafb8e
BASE="https://download.decode.is/s3/download?token=${TOKEN}&file="

# Exact basenames from the data-use bundle (note: deCODE/INTERMOUNTAIN/UKBB use _sumstat,
# FINNGEN + PDFF use _sumstats). HCC has no FinnGen arm.
FILES=(
  NAFL_deCODE_sumstat NAFL_INTERMOUNTAIN_sumstat NAFL_UKBB_sumstat NAFL_FINNGEN_sumstats
  Cirrhosis_deCODE_sumstat Cirrhosis_INTERMOUNTAIN_sumstat Cirrhosis_UKBB_sumstat Cirrhosis_FINNGEN_sumstats
  HCC_deCODE_sumstat HCC_INTERMOUNTAIN_sumstat HCC_UKBB_sumstat
  Proton_density_fat_fraction_UKBB_sumstats
)
b=${FILES[$((SLURM_ARRAY_TASK_ID-1))]}
echo "[task $SLURM_ARRAY_TASK_ID] file=$b  $(date)"

# md5 sidecar first (tiny)
curl -fsSL --retry 5 --connect-timeout 30 -o "${b}.tmd5sum" "${BASE}${b}.tmd5sum" || echo "  (no md5 sidecar for $b)"
exp=$(awk 'NR==1{print $1}' "${b}.tmd5sum" 2>/dev/null)

# skip if already downloaded and checksum matches
if [ -f "${b}.txt" ] && [ -n "${exp:-}" ]; then
  got=$(md5sum "${b}.txt" | awk '{print $1}')
  if [ "$exp" = "$got" ]; then echo "  SKIP (md5 OK): ${b}.txt"; exit 0; fi
  echo "  present but md5 mismatch -> re-downloading ${b}.txt"
fi

echo "  downloading ${b}.txt ..."
curl -fSL --retry 5 --connect-timeout 30 -o "${b}.txt" "${BASE}${b}.txt" || { echo "  FAILED download $b (token expired?)"; exit 1; }
got=$(md5sum "${b}.txt" | awk '{print $1}')
sz=$(ls -la "${b}.txt" | awk '{print $5}')
if [ -n "${exp:-}" ]; then
  [ "$exp" = "$got" ] && echo "  OK md5 $b (size=$sz)" || { echo "  MISMATCH $b exp=$exp got=$got (size=$sz)"; exit 2; }
else
  echo "  downloaded $b (no md5 to verify; size=$sz)"
fi
echo "[task $SLURM_ARRAY_TASK_ID] done $(date)"
