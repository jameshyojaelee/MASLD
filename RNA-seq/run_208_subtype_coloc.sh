#!/bin/bash
#SBATCH --job-name=208_subtype_coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/208_subtype_coloc_%j.out
#SBATCH --error=RNA-seq/logs/208_subtype_coloc_%j.err

# Activate rnaseq environment BEFORE set -u (ADDR2LINE is unbound in activate-binutils)
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

echo "=== 208_subtype_coloc.R ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

cd "$MASLD_PROJECT_ROOT"
mkdir -p RNA-seq/logs

Rscript RNA-seq/208_subtype_coloc.R

echo ""
echo "End: $(date)"
