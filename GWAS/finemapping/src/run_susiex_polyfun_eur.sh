#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=cpu,io
#SBATCH --job-name=susiexpolyfun
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/susiexpolyfun_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/susiexpolyfun_%A_%a.err

# SuSiEX EUR=PolyFun, EAS=1kg_eas (analog of run_susiex_topld_eur.sh).
# Output → output/{gwas}/<chr_pos>/SuSiEX_polyfun/

set -o pipefail
export PYTHONNOUSERSITE=1

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate susiex
set -u

LOCUS_ROW=${SLURM_ARRAY_TASK_ID:-1}
export LD_PANEL=polyfun
export SUSIEX_RESULTS_SUFFIX="_polyfun"
# EAS arm pinned to 1kg_eas (PolyFun is EUR-only)
export EAS_LD_DIR="${FM_DIR}/data/ld_ref/1kg_eas"

echo "[susiexpolyfun] LOCUS_ROW=${LOCUS_ROW}  EUR=polyfun  EAS=1kg_eas  start $(date)"
python src/10_run_susiex.py --locus-row "${LOCUS_ROW}" --threads 4
echo "[susiexpolyfun] LOCUS_ROW=${LOCUS_ROW}  done  $(date)"
