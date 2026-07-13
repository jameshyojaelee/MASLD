#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=meSuSiE
#SBATCH --output=GWAS/finemapping/logs/mesusie_mvp_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/mesusie_mvp_%A_%a.err
#
# WITHIN-MVP multi-ancestry MESuSiE: EUR (PolyFun LD) + EAS + AFR + AMR (1KG LD),
# per MVP trait, over MVP shared loci (results/susiex_mvp/shared_loci.csv from 09b).
# Missing strata (e.g. MVP_Cirrhosis_EAS) and <MIN_SNPS arms are dropped; MESuSiE
# runs on the surviving >= 2 arms.
# Output -> results/mesusie_mvp/<trait>/<locus_id>_mesusie.rds
#
# Submit (N = row count of results/susiex_mvp/shared_loci.csv minus header):
#   sbatch --qos=interactive --array=1-<N> src/run_mesusie_mvp.sh

set -eo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

LOCUS_ROW=${SLURM_ARRAY_TASK_ID:-1}

export ANCESTRIES="EUR,EAS,AFR,AMR"

# EUR uses PolyFun LD; non-EUR use 1KG (Berisa-Pickrell block .bim/.ld).
export LD_PANEL=polyfun
unset UKBB_LD_DIR
export EAS_LD_DIR="${FM_DIR}/data/ld_ref/1kg_eas"
export AFR_LD_DIR="${FM_DIR}/data/ld_ref/1kg_afr"
export AMR_LD_DIR="${FM_DIR}/data/ld_ref/1kg_amr"

export MESUSIE_RESULTS_SUFFIX="_mvp"
export SHARED_LOCI_FILE="${FM_DIR}/results/susiex_mvp/shared_loci.csv"

echo "[mesusie_mvp] LOCUS_ROW=${LOCUS_ROW}  ANCESTRIES=${ANCESTRIES}"
echo "[mesusie_mvp]   EUR=polyfun  EAS=1kg_eas  AFR=1kg_afr  AMR=1kg_amr"
echo "[mesusie_mvp]   start $(date)"
Rscript src/10b_run_mesusie.R
echo "[mesusie_mvp]   done  $(date)"
