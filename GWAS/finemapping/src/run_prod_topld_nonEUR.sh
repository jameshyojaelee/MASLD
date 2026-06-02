#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --array=0-241
#SBATCH --job-name=prod_topld_nonEUR
#SBATCH --output=GWAS/finemapping/logs/prod_topld_nonEUR_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/prod_topld_nonEUR_%A_%a.err

# 11 non-EUR GWAS × 22 chr = 242 SuSiE-COLOC tasks under LD_PANEL=topld.
# 5 EAS (2 from 2020 paper + 3 BBJ) + 3 PanUKBB AFR + 3 PanUKBB CSA (SAS).

set -o pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${BASE}"

GWAS_LIST=(
  2020_32514122_Cirrhosis_EAS   2020_32514122_HCC_EAS
  BBJ_ALT                        BBJ_AST                       BBJ_GGT
  PanUKBB_AFR_ALT                PanUKBB_AFR_AST               PanUKBB_AFR_GGT
  PanUKBB_CSA_ALT                PanUKBB_CSA_AST               PanUKBB_CSA_GGT
)

IDX=${SLURM_ARRAY_TASK_ID}
GWAS_IDX=$((IDX / 22))
CHR=$((IDX % 22 + 1))
GWAS=${GWAS_LIST[$GWAS_IDX]}

export LD_PANEL=topld
export COLOC_OUT_SUFFIX="_topld"
# Don't pin per-ancestry overrides; let get_ld_base_dir route to topld_<ancestry>.
unset EAS_LD_DIR AFR_LD_DIR SAS_LD_DIR

echo "[prod_topld_nonEUR] task ${IDX}: ${GWAS} chr${CHR}  start $(date)"
Rscript GWAS/finemapping/src/06_susie_coloc.R "${GWAS}" "${CHR}"
echo "[prod_topld_nonEUR] task ${IDX}: ${GWAS} chr${CHR}  done  $(date)"
