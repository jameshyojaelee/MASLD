#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=susiex
#SBATCH --output=logs/susiex_%j.out
#SBATCH --error=logs/susiex_%j.err
#
# Run SuSiEX joint EUR+EAS fine-mapping for one locus.
#
# Usage (single locus by row number, 1-indexed):
#   LOCUS_ROW=1 sbatch 10_run_susiex.sh
#
# Usage (SLURM array over all loci):
#   N=$(tail -n +2 results/susiex/shared_loci.csv | wc -l)
#   sbatch --array=1-${N} 10_run_susiex.sh
#
# Prerequisites: run 09_identify_shared_loci.R first to generate shared_loci.csv.
# Environment: susiex conda env (python 3.10 + pandas/numpy/scipy)
# SuSiEX binary: bin/SuSiEx (static, copied from repo bin_static/)

set -euo pipefail

# Prevent user-level Python packages from shadowing conda env packages
export PYTHONNOUSERSITE=1

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

# Activate conda environment
eval "$(micromamba shell hook -s bash)"
micromamba activate susiex

# Determine locus row from SLURM array task ID or LOCUS_ROW env var
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    LOCUS_ROW=${SLURM_ARRAY_TASK_ID}
else
    LOCUS_ROW=${LOCUS_ROW:-1}
fi

echo "Running SuSiEX for locus row: ${LOCUS_ROW}"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Start: $(date)"

python src/10_run_susiex.py --locus-row "${LOCUS_ROW}" --threads 4

echo "Done: $(date)"
