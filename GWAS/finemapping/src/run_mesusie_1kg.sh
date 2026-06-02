#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=mesusie1kg
#SBATCH --output=GWAS/finemapping/logs/mesusie1kg_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/mesusie1kg_%A_%a.err
#
# Run MESuSiE multi-ancestry fine-mapping with LD_PANEL=1kg for the same 140
# shared EUR×EAS loci. Pure 1KG baseline (no sghatan dependency).
#
# Usage (SLURM array over all shared loci):
#   sbatch --array=1-140 src/run_mesusie_1kg.sh

set -eo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

LOCUS_ROW=${SLURM_ARRAY_TASK_ID:-1}
export LD_PANEL=1kg
export MESUSIE_RESULTS_SUFFIX="_1kg"
# Leave EAS_LD_DIR / EUR-side overrides UNSET so dispatch resolves to 1kg_<ancestry>
unset EAS_LD_DIR AFR_LD_DIR SAS_LD_DIR UKBB_LD_DIR

echo "[mesusie1kg] LOCUS_ROW=${LOCUS_ROW}  LD_PANEL=${LD_PANEL}  start $(date)"
Rscript src/10b_run_mesusie.R
echo "[mesusie1kg] LOCUS_ROW=${LOCUS_ROW}  done  $(date)"
