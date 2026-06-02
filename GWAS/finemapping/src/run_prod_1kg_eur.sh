#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --array=0-373
#SBATCH --job-name=prod_1kg_eur
#SBATCH --output=GWAS/finemapping/logs/prod_1kg_eur_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/prod_1kg_eur_%A_%a.err

# Production re-run: 17 EUR GWAS × 22 chromosomes = 374 SuSiE-COLOC tasks
# using LD_PANEL=1kg. Output goes to results/susie_coloc_1kg/<gwas>/ via
# COLOC_OUT_SUFFIX env var (v1 UKBB snapshot in results/susie_coloc_ukbb_sghatan_v1_2026-04-21/
# is preserved; current results/susie_coloc/ is also preserved).

set -o pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${BASE}"

# Decode array index: 17 GWAS × 22 chr
GWAS_LIST=(
  2019_31311600_NAFLD_EUR
  2020_32298765_NAFLD_EUR
  2021_34128465_PDFF_EUR
  2021_34841290_NAFLD_EUR
  2021_34957434_PDFF_EUR
  2022_36402844_PDFF_EUR
  2023_36280732_NAFLD_deCode_EUR
  2023_36280732_NAFLD_Intermountain_EUR
  2023_36280732_NAFLD_UKBB_EUR
  UKBB_ALT
  UKBB_AST
  UKBB_GGT
  FinnGen_NAFLD
  FinnGen_NASH
  FinnGen_HCC
  Ghouse_Cirrhosis
  Ghouse_HCC
)

IDX=${SLURM_ARRAY_TASK_ID}
GWAS_IDX=$((IDX / 22))
CHR=$((IDX % 22 + 1))
GWAS=${GWAS_LIST[$GWAS_IDX]}

export LD_PANEL=1kg
export COLOC_OUT_SUFFIX="_1kg"

echo "[prod_1kg] task ${IDX}: ${GWAS} chr${CHR}  start $(date)"
Rscript GWAS/finemapping/src/06_susie_coloc.R "${GWAS}" "${CHR}"
echo "[prod_1kg] task ${IDX}: ${GWAS} chr${CHR}  done  $(date)"
