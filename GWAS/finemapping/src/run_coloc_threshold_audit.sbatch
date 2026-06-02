#!/bin/bash
#SBATCH --job-name=coloc_threshold_audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/coloc_threshold_audit_%j.out
#SBATCH --error=logs/coloc_threshold_audit_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

echo "=== COLOC Threshold Audit ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript src/diagnose_coloc_thresholds.R 2>&1 | tee logs/coloc_threshold_audit_${SLURM_JOB_ID}.log

echo ""
echo "Done: $(date)"
