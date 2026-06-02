#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=mesusietopld
#SBATCH --output=GWAS/finemapping/logs/mesusietopld_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/mesusietopld_%A_%a.err
#
# Run MESuSiE with EUR arm = topld_eur, EAS arm = 1kg_eas (mirrors the SuSiEX
# "topld" run — EUR-only LD swap). Output → results/mesusie_topld/.
#
# Usage (SLURM array over all shared loci):
#   sbatch --array=1-140 src/run_mesusie_topld.sh

set -eo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

LOCUS_ROW=${SLURM_ARRAY_TASK_ID:-1}
export LD_PANEL=topld
export MESUSIE_RESULTS_SUFFIX="_topld"
# EAS arm should NOT swap to topld_eas in this variant — pin to 1kg_eas.
export EAS_LD_DIR="${FM_DIR}/data/ld_ref/1kg_eas"

echo "[mesusietopld] LOCUS_ROW=${LOCUS_ROW}  LD_PANEL=${LD_PANEL}  EAS_LD_DIR=${EAS_LD_DIR}  start $(date)"
Rscript src/10b_run_mesusie.R
echo "[mesusietopld] LOCUS_ROW=${LOCUS_ROW}  done  $(date)"
