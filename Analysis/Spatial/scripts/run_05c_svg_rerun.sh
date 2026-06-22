#!/bin/bash
#SBATCH --job-name=squidpy
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/05c_svg_rerun_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/05c_svg_rerun_%j.err
#
# M8 remediation re-run: 05c_spatially_variable_genes.py was edited to
# DONOR_COL='individual' (collapses GSE192741 H35=JBO014+JBO015 to one biological
# replicate -> Steatotic n=2 not 3) but had NOT been re-run. The cited SVG counts
# (Healthy 370 / Steatotic 365 / disease_emergent 94) come from the pre-fix
# sample_id build (results/svg/, mtime 2026-06-01). This re-run refreshes
# svgs_Healthy.csv, svgs_Steatotic.csv, differential_svgs.csv under the corrected
# individual-level donor blocking. Input h5ad obs verified to carry the
# 'individual' column (H35/H36/H37/H38).
set -euo pipefail

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ"

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-NA}"
micromamba run -n spatial python "$PROJ/Analysis/Spatial/scripts/05c_spatially_variable_genes.py"
echo "[$(date)] 05c SVG re-run complete"
