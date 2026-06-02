#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --partition=cpu,io,dev
#SBATCH --array=0-373
#SBATCH --job-name=polyfun-susie-coloc
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/prod_polyfun_eur_portable_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/prod_polyfun_eur_portable_%A_%a.err

# Phase 9d portable variant — runnable by any user with read access to
# /gpfs/commons/home/jameslee. Uses absolute paths to micromamba so
# nsaravanan / tzhong / etc. don't need micromamba initialized in their
# shell. Idempotent skip handles overlap with jameslee's parallel queue.

set -o pipefail

export MAMBA_EXE="/gpfs/commons/home/jameslee/.local/bin/micromamba"
export MAMBA_ROOT_PREFIX="/gpfs/commons/home/jameslee/micromamba"
eval "$(${MAMBA_EXE} shell hook -s bash)"
micromamba activate rnaseq

if [ -z "${CONDA_PREFIX:-}" ]; then
  echo "ERROR: micromamba activation failed. Check that ${MAMBA_EXE} is readable."
  exit 1
fi

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

echo "[prod_polyfun_portable] task=${IDX} user=$(whoami) gwas=${GWAS} chr=${CHR} start=$(date)"

if [ -f "${OUT_CSV}" ]; then
  SIZE=$(stat -c %s "${OUT_CSV}")
  if [ "${SIZE}" -gt 1000 ]; then
    echo "[prod_polyfun_portable] task=${IDX} output exists (${SIZE} bytes) — skipping"
    exit 0
  fi
fi

Rscript GWAS/finemapping/src/06_susie_coloc.R "${GWAS}" "${CHR}"
rc=$?
echo "[prod_polyfun_portable] task=${IDX} ${GWAS} chr${CHR} done=$(date) rc=${rc}"
exit ${rc}
