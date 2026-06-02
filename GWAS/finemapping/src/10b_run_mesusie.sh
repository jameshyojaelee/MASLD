#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=mesusie
#SBATCH --output=logs/mesusie_%j.out
#SBATCH --error=logs/mesusie_%j.err
#
# Run MESuSiE multi-ancestry fine-mapping for one locus.
#
# Usage (SLURM array over all shared loci):
#   N=$(tail -n +2 results/susiex/shared_loci.csv | wc -l)
#   sbatch --array=1-${N} src/10b_run_mesusie.sh

set -eo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "Running MESuSiE for locus row: ${SLURM_ARRAY_TASK_ID:-${LOCUS_ROW:-1}}"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Start: $(date)"

Rscript src/10b_run_mesusie.R

echo "Done: $(date)"
