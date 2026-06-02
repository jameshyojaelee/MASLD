#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --array=0-373
#SBATCH --job-name=prod_topld_eur
#SBATCH --output=GWAS/finemapping/logs/prod_topld_eur_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/prod_topld_eur_%A_%a.err

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

export LD_PANEL=topld
export COLOC_OUT_SUFFIX="_topld"

echo "[prod_topld] task ${IDX}: ${GWAS} chr${CHR}  start $(date)"
Rscript GWAS/finemapping/src/06_susie_coloc.R "${GWAS}" "${CHR}"
echo "[prod_topld] task ${IDX}: ${GWAS} chr${CHR}  done  $(date)"
