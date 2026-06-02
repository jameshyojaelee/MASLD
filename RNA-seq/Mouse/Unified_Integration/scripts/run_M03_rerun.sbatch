#!/bin/bash
#SBATCH --job-name=metafor
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/M03_rerun_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/M03_rerun_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --qos=nslab

# M03 re-run: meta-analysis + pooled dream with two bug fixes (2026-05-26):
#  1. eBayes(fit) removed after dream() — was double variance shrinkage
#  2. Meta-analysis SE now reads SE_unmoderated from M02 output instead of logFC/t

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration

echo "=== M03 Re-run: Meta-Analysis + Pooled Dream (bug fixes) ==="
echo "Start: $(date)"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo ""

Rscript scripts/M03_mouse_meta_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: M03"; exit 1; fi

echo ""
echo "=== M03 complete: $(date) ==="
