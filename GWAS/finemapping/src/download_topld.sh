#!/bin/bash -l
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --partition=cpu
#SBATCH --array=0-3
#SBATCH --job-name=download_topld
#SBATCH --output=logs/download_topld_%A_%a.out
#SBATCH --error=logs/download_topld_%A_%a.err

# download_topld.sh
# Download per-chromosome TOP-LD R²/D' files from
# http://topld.genetics.unc.edu/downloads/topld_v1/{POP}/SNV/
# One array task per ancestry (EUR/AFR/EAS/SAS).
#
# File pattern: {POP}_chr{N}_Dprime_extended_0.2_1000000.csv.gz
#               {POP}_chr{N}_Dprime_extended_0.2_1000000_info.csv.gz
# Threshold note: "0.2" = R² ≥ 0.2 cutoff; "1000000" = 1Mb max distance window.
# Total raw download per ancestry: ~7-15 GB compressed. ~88 files total.

set -euo pipefail

POPS=(EUR AFR EAS SAS)
POP=${POPS[$SLURM_ARRAY_TASK_ID]}

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT_DIR="${BASE}/GWAS/finemapping/data/ld_ref/topld_raw/${POP}"
URL_BASE="http://topld.genetics.unc.edu/downloads/topld_v1/${POP}/SNV"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

echo "============================================================"
echo "Downloading TOP-LD ${POP} (chr1-22) starting $(date)"
echo "============================================================"

N_OK=0
N_SKIP=0
N_FAIL=0

for chr in {1..22}; do
  for suffix in "" "_info"; do
    fname="${POP}_chr${chr}_Dprime_extended_0.2_1000000${suffix}.csv.gz"
    url="${URL_BASE}/${fname}"

    # Skip if already downloaded and non-empty
    if [ -s "${fname}" ]; then
      N_SKIP=$((N_SKIP + 1))
      continue
    fi

    # wget -c for resume support
    if wget -q --tries=3 --timeout=60 -c "${url}" -O "${fname}.tmp"; then
      mv "${fname}.tmp" "${fname}"
      sz=$(stat -c%s "${fname}")
      sz_mb=$((sz / 1024 / 1024))
      echo "  [${POP}] chr${chr}${suffix}: ${sz_mb} MB ✓"
      N_OK=$((N_OK + 1))
    else
      rm -f "${fname}.tmp"
      echo "  [${POP}] chr${chr}${suffix}: FAILED"
      N_FAIL=$((N_FAIL + 1))
    fi
  done
done

echo ""
echo "[${POP}] Done $(date): ${N_OK} downloaded, ${N_SKIP} skipped (pre-existing), ${N_FAIL} failed"
echo "Total disk: $(du -sh "${OUT_DIR}" | cut -f1)"
