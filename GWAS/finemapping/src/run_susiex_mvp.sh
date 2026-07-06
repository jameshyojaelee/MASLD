#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu,io
#SBATCH --job-name=SuSiEx
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/susiex_mvp_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/susiex_mvp_%A_%a.err
#
# WITHIN-MVP multi-ancestry SuSiEX: EUR + EAS + AFR + AMR (all on 1KG LD), per
# MVP trait, over MVP shared loci (results/susiex_mvp/shared_loci.csv from 09b).
# Per-locus arms with a missing stratum (e.g. MVP_Cirrhosis_EAS) or < MIN_SNPS
# are dropped automatically; SuSiEX runs on the surviving >= 2 arms.
# Output -> results/susiex_mvp/<trait>/<locus_id>.{summary,cs,snp}
#
# Submit (N = row count of results/susiex_mvp/shared_loci.csv minus header):
#   sbatch --qos=interactive --array=1-<N> src/run_susiex_mvp.sh

set -o pipefail
export PYTHONNOUSERSITE=1

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate susiex
set -u

LOCUS_ROW=${SLURM_ARRAY_TASK_ID:-1}

# MVP ancestries (EUR anchor first). SAS is not an MVP stratum.
export ANCESTRIES="EUR,EAS,AFR,AMR"

# All arms on 1KG (SuSiEX needs per-ancestry genotypes; PolyFun ships LD only).
export LD_PANEL=1kg
unset UKBB_LD_DIR EAS_LD_DIR AFR_LD_DIR AMR_LD_DIR SAS_LD_DIR

export SUSIEX_RESULTS_SUFFIX="_mvp"
export SHARED_LOCI_FILE="${FM_DIR}/results/susiex_mvp/shared_loci.csv"

echo "[susiex_mvp] LOCUS_ROW=${LOCUS_ROW}  ANCESTRIES=${ANCESTRIES}"
echo "[susiex_mvp]   EUR=1kg_eur  EAS=1kg_eas  AFR=1kg_afr  AMR=1kg_amr"
echo "[susiex_mvp]   start $(date)"
python src/10_run_susiex.py --locus-row "${LOCUS_ROW}" --threads 4
echo "[susiex_mvp]   done  $(date)"
