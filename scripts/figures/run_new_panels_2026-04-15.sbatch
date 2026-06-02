#!/bin/bash
#SBATCH --job-name=fig_new_panels
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig_new_panels_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig_new_panels_%j.err

set -eo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE"
mkdir -p "$BASE/scripts/figures/logs"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

echo "[$(date)] === Fig 1: % patients vs Log2FC cutoff heatmap ==="
Rscript "$BASE/scripts/figures/fig1_patient_lfc_cutoff.R"

echo "[$(date)] === Fig 2: MASLD characterization bar ==="
Rscript "$BASE/scripts/figures/fig2_masld_characterization.R"

echo "[$(date)] === Fig 3: Cross-ancestry COLOC elastic net classifier ==="
Rscript "$BASE/scripts/figures/fig3_cross_ancestry_classifier.R"

echo "[$(date)] === All 3 new panels complete ==="
