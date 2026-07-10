#!/bin/bash
# Render the up-to-date per-cell-type Hotspot autocorrelation module heatmaps
# into a Fig 3 panels subdir, restricted to the 5 protocol-clean cell types
# (global + endothelial_cells were deferred at the 2026-05-22 remediation).
#SBATCH --job-name=hotspot
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/hotspot_modules/logs/figS_heatmap_fig3_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/hotspot_modules/logs/figS_heatmap_fig3_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

export HOTSPOT_FIG_OUT="figures/main/fig3_RNAseq/panels/hotspot_module_heatmaps"
export HOTSPOT_RUN_CTS="hepatocytes,macrophages,fibroblasts,cholangiocytes,tcells"

micromamba run -n rapids_singlecell python scripts/figures/figS_hotspot_heatmap.py
