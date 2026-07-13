#!/bin/bash
#SBATCH --job-name=plotgardener
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=500G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --output=logs/plotgardener_%j.out
#SBATCH --error=logs/plotgardener_%j.err
# Whole-library Cas13 guide-track coverage figure.
#   1) build_all_guide_coords.py  -> data/guide_qc_all_library.tsv   (every targeting
#      guide in data/guides/library.csv; non-targeting controls excluded)
#   2) spotcheck_guides_plotgardener.R -> figures/guide_qc/spotcheck2_gene_tracks_plotgardener.pdf
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Cas13_Library_Design
mkdir -p logs
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "[1/2] building all-library coordinate table ..."
python scripts/guides/build_all_guide_coords.py

echo "[2/2] rendering plotgardener whole-library tracks ..."
GENES_PER_PAGE="${GENES_PER_PAGE:-6}" Rscript scripts/guides/spotcheck_guides_plotgardener.R

echo "DONE"
