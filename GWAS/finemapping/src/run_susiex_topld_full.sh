#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=susiextopldfull
#SBATCH --output=GWAS/finemapping/logs/susiextopldfull_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/susiextopldfull_%A_%a.err

# Re-run SuSiEX joint EUR+EAS fine-mapping with BOTH arms on TOP-LD (parallel
# to run_susiex_topld_eur.sh which kept the EAS arm on 1kg_eas). Closes the
# symmetry across the LD-panel comparison axes.
#
# Prerequisite: chr-level PLINK files for topld_eas built via
# build_topld_eas_chrlevel.sh (chr<N>_eas.{bed,bim,fam}). SuSiEX EAS arm reads
# chr-level PLINK (10_run_susiex.py:get_eas_ref_prefix lines 352-356).

set -o pipefail
export PYTHONNOUSERSITE=1

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate susiex
set -u

LOCUS_ROW=${SLURM_ARRAY_TASK_ID:-1}
export LD_PANEL=topld
export SUSIEX_RESULTS_SUFFIX="_topld_full"
# EAS arm SWAPS to topld_eas (the 7c new build). Set explicitly to avoid any
# fallback-to-1kg_eas behavior in the dispatcher.
export EAS_LD_DIR="${FM_DIR}/data/ld_ref/topld_eas"

echo "[susiextopldfull] LOCUS_ROW=${LOCUS_ROW}  LD_PANEL=${LD_PANEL}  EAS_LD_DIR=${EAS_LD_DIR}  start $(date)"
python src/10_run_susiex.py --locus-row "${LOCUS_ROW}" --threads 4
echo "[susiextopldfull] LOCUS_ROW=${LOCUS_ROW}  done  $(date)"
