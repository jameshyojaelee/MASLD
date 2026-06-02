#!/bin/bash
#SBATCH --job-name=limma-voom
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/M02_rerun_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/M02_rerun_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --qos=nslab

# M02 re-run: per-diet DE with SE_unmoderated + df_total columns added
# Bug fix 2026-05-26: M02 now saves unmoderated SE (stdev.unscaled * sigma)
# so M03 meta-analysis uses proper SE instead of moderated logFC/t.

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration

echo "=== M02 Re-run: Per-Diet DE (SE fix) ==="
echo "Start: $(date)"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo ""

Rscript scripts/M02_mouse_per_diet_de.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: M02"; exit 1; fi

echo ""
echo "=== M02 complete: $(date) ==="
