#!/bin/bash
#SBATCH --job-name=hotspot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/hotspot_modules/logs/511_donorcollapse_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/hotspot_modules/logs/511_donorcollapse_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs

echo "=========================================================="
echo "[1/3] pooled-cell (weighted) donor scores from parquet"
echo "=========================================================="
micromamba run -n rnaseq python \
  Analysis/SingleCell/scripts/hotspot_modules/511_pooled_cell_donor_scores.py

echo "=========================================================="
echo "[2/3] donor-level disease-stage recompute + all_modules_donor.tsv"
echo "=========================================================="
micromamba run -n rnaseq Rscript \
  Analysis/SingleCell/scripts/hotspot_modules/511_donorlevel_disease_stage.R

echo "=========================================================="
echo "[3/3] regenerate Fig3G module heatmap (donor-level, NEW filenames)"
echo "=========================================================="
micromamba run -n rnaseq Rscript \
  scripts/figures/singlecell_module_heatmap_donorlevel.R

echo "[DONE] donor-collapse hotspot recompute complete"
