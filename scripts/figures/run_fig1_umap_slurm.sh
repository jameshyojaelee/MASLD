#!/bin/bash
#SBATCH --job-name=liver_fig1_umap
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/fig1_umap_%j.log

# Re-compute Human + Mouse UMAPs, then regenerate Fig 1
set -euo pipefail

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$MASLD_PROJECT_ROOT"
mkdir -p logs figures

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== Step 1: Human UMAP (10 cohorts) ==="
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/09_batch_correction_umap.R

echo "=== Step 2: Mouse UMAP (7 datasets) ==="
Rscript RNA-seq/Mouse/Unified_Integration/scripts/M06_mouse_batch_correction_umap.R

echo "=== Step 3: Generate Fig 1 ==="
Rscript scripts/figures/fig1.R

echo "=== Done ==="
ls -la figures/fig1.pdf
