#!/bin/bash -l
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=pf11wk
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/pf11wk_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/pf11wk_%A_%a.err
#
# 32-worker parallel SuSiE-COLOC for chr11 / 2022_36402844_PDFF_EUR (PolyFun LD).
# Workers stride: each handles eGenes where (i-1) %% 32 == WORKER_ID.
# All workers read shared checkpoint_chr11.csv (1001 done genes) and skip those.
# Each writes to worker_chr11_w<id>.csv; merge_polyfun_chr11_workers.R combines them.
#
# Memory: --mem=24G gives ~70% headroom over observed peak (~14 GB) per the
# original packed wrapper's note. Safe for chr11's densest LD blocks.
#
# Submit:
#   sbatch --array=0-31 src/run_polyfun_chr11_workers.sh

set -o pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${BASE}"

export LD_PANEL=polyfun
export COLOC_OUT_SUFFIX="_polyfun"
unset UKBB_LD_DIR EAS_LD_DIR AFR_LD_DIR SAS_LD_DIR

export WORKER_ID=${SLURM_ARRAY_TASK_ID:-0}
export N_WORKERS=32

GWAS=2022_36402844_PDFF_EUR
CHR=11

echo "[pf11wk] worker ${WORKER_ID}/${N_WORKERS}  GWAS=${GWAS} chr${CHR}  start $(date)"
echo "[pf11wk] node $(hostname)  cpus=${SLURM_CPUS_ON_NODE:-?}"

Rscript GWAS/finemapping/src/06_susie_coloc.R "${GWAS}" "${CHR}" 2>&1 \
  | sed "s/^/[w${WORKER_ID}] /"

rc=$?
echo "[pf11wk] worker ${WORKER_ID} done $(date) rc=${rc}"
exit $rc
