#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=mesusie4way
#SBATCH --output=GWAS/finemapping/logs/mesusie4way_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/mesusie4way_%A_%a.err
#
# 4-ancestry MESuSiE: EUR (PolyFun) + EAS (1kg_eas) + AFR (1kg_afr) + SAS (1kg_sas)
# for 140 shared loci × ALT/AST/GGT.
# Output → results/mesusie_4way/<trait>/<locus_id>_mesusie.rds
#
# Submit:
#   sbatch --array=1-140 src/run_mesusie_4way.sh

set -eo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

LOCUS_ROW=${SLURM_ARRAY_TASK_ID:-1}

export ANCESTRIES="EUR,EAS,AFR,SAS"

export LD_PANEL=polyfun
unset UKBB_LD_DIR
export EAS_LD_DIR="${FM_DIR}/data/ld_ref/1kg_eas"
export AFR_LD_DIR="${FM_DIR}/data/ld_ref/1kg_afr"
export SAS_LD_DIR="${FM_DIR}/data/ld_ref/1kg_sas"

export MESUSIE_RESULTS_SUFFIX="_4way"

echo "[mesusie4way] LOCUS_ROW=${LOCUS_ROW}  ANCESTRIES=${ANCESTRIES}"
echo "[mesusie4way]   EUR=polyfun  EAS=1kg_eas  AFR=1kg_afr  SAS=1kg_sas"
echo "[mesusie4way]   start $(date)"
Rscript src/10b_run_mesusie.R
echo "[mesusie4way]   done  $(date)"
