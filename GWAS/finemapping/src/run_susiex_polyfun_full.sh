#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=cpu,io
#SBATCH --job-name=susiexpolyfunfull
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/susiexpolyfunfull_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/susiexpolyfunfull_%A_%a.err

# SuSiEX EUR=PolyFun, EAS=topld_eas (closes symmetry with TOP-LD-full).
# Output → output/{gwas}/<chr_pos>/SuSiEX_polyfun_full/

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
export SUSIEX_RESULTS_SUFFIX="_polyfun_full"
# EAS arm uses TOP-LD-EAS (the chr-level PLINK we built in 7c)
export EAS_LD_DIR="${FM_DIR}/data/ld_ref/topld_eas"

echo "[susiexpolyfunfull] LOCUS_ROW=${LOCUS_ROW}  EUR=polyfun  EAS=topld_eas  start $(date)"
python src/10_run_susiex.py --locus-row "${LOCUS_ROW}" --threads 4
echo "[susiexpolyfunfull] LOCUS_ROW=${LOCUS_ROW}  done  $(date)"
