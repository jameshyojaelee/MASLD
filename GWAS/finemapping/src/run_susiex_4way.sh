#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=cpu,io
#SBATCH --job-name=susiex4way
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/susiex4way_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/susiex4way_%A_%a.err
#
# 4-ancestry SuSiEX: EUR (PolyFun) + EAS (1kg_eas) + AFR (1kg_afr) + SAS (1kg_sas)
# for 140 shared loci × ALT/AST/GGT.
# Output → results/susiex_4way/<trait>/<locus_id>.{summary,cs,snp}
#
# Submit:
#   sbatch --array=1-140 src/run_susiex_4way.sh

set -o pipefail
export PYTHONNOUSERSITE=1

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate susiex
set -u

LOCUS_ROW=${SLURM_ARRAY_TASK_ID:-1}

# Ancestry list (order = SuSiEX arm order; EUR first as anchor)
export ANCESTRIES="EUR,EAS,AFR,SAS"

# All ancestries on 1KG (SuSiEX needs per-ancestry .bed/.fam genotypes; PolyFun
# ships LD matrices only, so PolyFun EUR cannot be used for SuSiEX).
export LD_PANEL=1kg
unset UKBB_LD_DIR EAS_LD_DIR AFR_LD_DIR SAS_LD_DIR

export SUSIEX_RESULTS_SUFFIX="_4way"

echo "[susiex4way] LOCUS_ROW=${LOCUS_ROW}  ANCESTRIES=${ANCESTRIES}"
echo "[susiex4way]   EUR=1kg_eur  EAS=1kg_eas  AFR=1kg_afr  SAS=1kg_sas"
echo "[susiex4way]   start $(date)"
python src/10_run_susiex.py --locus-row "${LOCUS_ROW}" --threads 4
echo "[susiex4way]   done  $(date)"
