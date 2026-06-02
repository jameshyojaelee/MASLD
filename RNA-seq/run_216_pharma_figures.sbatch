#!/bin/bash
#SBATCH --job-name=s216_pharma_fig
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/216_pharma_figures_%j.out
#SBATCH --error=logs/216_pharma_figures_%j.err

# Environment BEFORE set -u (activate-binutils references unbound ADDR2LINE)
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

echo "=== Script 216: Pharmacogenomic Synthesis Figures ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/RNA-seq"

mkdir -p logs

Rscript 216_pharma_figures.R

echo "End: $(date)"
echo "=== Script 216 complete ==="
