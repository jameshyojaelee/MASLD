#!/bin/bash
#SBATCH --job-name=liver_mouse_viz
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/viz_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/viz_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu

set -euo pipefail

SCRIPTS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/scripts"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "===== M05: Volcano Plots ====="
Rscript "${SCRIPTS}/M05_mouse_volcano_plots.R"

echo ""
echo "===== M06: Batch Correction + UMAP ====="
Rscript "${SCRIPTS}/M06_mouse_batch_correction_umap.R"

echo ""
echo "===== M07: Pathway Analysis (GSEA) ====="
Rscript "${SCRIPTS}/M07_mouse_pathway_analysis.R"

echo ""
echo "===== All mouse visualization scripts complete ====="
