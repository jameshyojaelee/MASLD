#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --partition=cpu,io,dev,bigmem
#SBATCH --array=0-373
#SBATCH --job-name=polyfun-susie-coloc
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/prod_polyfun_eur_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/prod_polyfun_eur_%A_%a.err

# Phase 9d — 17 EUR GWAS × 22 chr = 374 SuSiE-COLOC tasks under
# LD_PANEL=polyfun (rebuilt by Phase 9 build_polyfun_blocks.py with
# single-window-per-block LD).
#
# Multi-partition cpu+io+dev+bigmem for maximum scheduling spread.
# 16G memory request: 8d MaxRSS history showed peak ~14G; 16G adds
# safety margin without bloating allocation footprint.

set -o pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${BASE}"

GWAS_LIST=(
  2019_31311600_NAFLD_EUR  2020_32298765_NAFLD_EUR  2021_34128465_PDFF_EUR
  2021_34841290_NAFLD_EUR  2021_34957434_PDFF_EUR  2022_36402844_PDFF_EUR
  2023_36280732_NAFLD_deCode_EUR  2023_36280732_NAFLD_Intermountain_EUR
  2023_36280732_NAFLD_UKBB_EUR  UKBB_ALT  UKBB_AST  UKBB_GGT
  FinnGen_NAFLD  FinnGen_NASH  FinnGen_HCC  Ghouse_Cirrhosis  Ghouse_HCC
)

IDX=${SLURM_ARRAY_TASK_ID}
GWAS_IDX=$((IDX / 22))
CHR=$((IDX % 22 + 1))
GWAS=${GWAS_LIST[$GWAS_IDX]}

export LD_PANEL=polyfun
export COLOC_OUT_SUFFIX="_polyfun"
unset UKBB_LD_DIR EAS_LD_DIR AFR_LD_DIR SAS_LD_DIR

OUT_CSV="${BASE}/GWAS/finemapping/results/susie_coloc_polyfun/${GWAS}/susie_coloc_chr${CHR}.csv"

echo "[prod_polyfun] task ${IDX}: ${GWAS} chr${CHR}  start $(date)"

# Idempotency: skip if CART (or prior run) already produced this output
if [ -f "${OUT_CSV}" ]; then
  SIZE=$(stat -c %s "${OUT_CSV}")
  if [ "${SIZE}" -gt 1000 ]; then
    echo "[prod_polyfun] task ${IDX}: output exists (${SIZE} bytes) — skipping"
    exit 0
  fi
fi

Rscript GWAS/finemapping/src/06_susie_coloc.R "${GWAS}" "${CHR}"
echo "[prod_polyfun] task ${IDX}: ${GWAS} chr${CHR}  done  $(date)"
