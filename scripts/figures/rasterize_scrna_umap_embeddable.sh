#!/bin/bash
#SBATCH --job-name=rasterize_scrna_umap
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --output=scripts/figures/logs/rasterize_scrna_umap_%j.log

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PANEL_DIR="${BASE}/figures/main/fig3_RNAseq/panels"
SRC="${PANEL_DIR}/scrna_umap.pdf"
OUT="${PANEL_DIR}/fig3f_scrna_umap_embeddable.pdf"
TMP=$(mktemp -d)

echo "[$(date)] rasterizing ${SRC} → ${OUT}"
echo "[tmp] ${TMP}"

# Step 1: rasterize each page to PNG at 300 DPI via pdftoppm
pdftoppm -r 300 -png "${SRC}" "${TMP}/page"

echo "[$(date)] rasterized pages:"
ls -la "${TMP}/page"*.png

# Step 2: wrap PNG(s) back into a PDF using ghostscript
# This produces a fully rasterized PDF with no vector content
PAGE_COUNT=$(ls "${TMP}/page"*.png | wc -l)
echo "[page count] ${PAGE_COUNT}"

gs -dNOPAUSE -dBATCH -dQUIET \
   -sDEVICE=pdfwrite \
   -dPDFSETTINGS=/prepress \
   -dCompatibilityLevel=1.4 \
   -dAutoRotatePages=/None \
   -sOutputFile="${OUT}" \
   "${TMP}/page"*.png

echo "[$(date)] done → ${OUT} ($(du -sh "${OUT}" | cut -f1))"
ls -la "${OUT}"

rm -rf "${TMP}"
