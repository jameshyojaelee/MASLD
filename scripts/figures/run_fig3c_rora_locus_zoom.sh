#!/bin/bash
#SBATCH --job-name=fig3c_rora
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig3c_rora_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig3c_rora_%j.err
#
# fig3c: GGT chr15 / RORA locus LD-zoom panel (the "external" c in the Fig 3 composite).
# Single-step pipeline — figS09_locus_zoom_pg.R (plotgardener rewrite) writes
# the LD-zoom plot DIRECTLY to
# fig3_regulatory_architecture/panels/rora_locus_zoom_ggt_chr15.pdf
# via the 4th CLI arg (output override). No copy step needed.
#
# 2026-04-29 — switched from rnaseq+figS09_locus_zoom.R (ggplot/patchwork) to
# plotgardener+figS09_locus_zoom_pg.R for coordinate-perfect tracks (native
# plotGenes / plotGenomeLabel) and prettier genomic axis. The legacy ggplot
# script is preserved at scripts/figures/figS09_locus_zoom.R and remains the
# producer of figS09 panels.

eval "$(micromamba shell hook --shell bash)"
micromamba activate plotgardener
set -uo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "${BASE}/scripts/figures/logs"
cd "${BASE}"

echo "[fig3c] start $(date) — RORA locus (GGT chr15:60883281), plotgardener path"
Rscript scripts/figures/figS09_locus_zoom_pg.R \
  locus_GGT_chr15_60883281 \
  GGT \
  RORA \
  "${BASE}/figures/main/fig3_regulatory_architecture/panels/rora_locus_zoom_ggt_chr15.pdf"
echo "[fig3c] done  $(date)"
