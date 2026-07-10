#!/bin/bash
#SBATCH --job-name=plotgardener
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/locuszoom_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/locuszoom_%j.err
#
# Consolidated Fig2 locus-zoom regen (2026-07-01). Renders the full plotgardener
# roster into the single figures/main/fig2_genetics/locus_zoom/ dir with the
# harmonized <TRAIT>_<GENE>.pdf naming and the compacted (~5.4x4.2-5.1in, uniform
# 6pt font floor) layout. Absorbs the former run_locuszoom_secondary.sh and
# run_fig3c_rora_locus_zoom.sh (both retired — fully covered here).
#
#   - enzyme / cross-ancestry panels (figS09_locus_zoom_pg.R):
#     title "<gene> · <TRAIT> (cross-ancestry)" + GWAS tracks "<TRAIT> · EUR/EAS"
#   - disease / PDFF panels (figS09_locus_zoom_disease.R): GWAS track "<TRAIT> · EUR"
#
# CYFIP2/CYP26A1 (ggplot-based, different renderer) regen separately via
# run_locuszoom_cyfip2_cyp26a1.sh (rnaseq env).
set +e   # one render failing must not abort the rest
eval "$(micromamba shell hook --shell bash)"
micromamba activate plotgardener
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE"
P="$BASE/figures/main/fig2_genetics/locus_zoom"
mkdir -p "$P" scripts/figures/logs

echo "===== ENZYME (cross-ancestry, GGT) ====="
Rscript scripts/figures/figS09_locus_zoom_pg.R locus_GGT_chr15_60883281  GGT RORA   "$P/GGT_RORA.pdf"
Rscript scripts/figures/figS09_locus_zoom_pg.R locus_GGT_chr1_16505320   GGT EPHA2  "$P/GGT_EPHA2.pdf"
Rscript scripts/figures/figS09_locus_zoom_pg.R locus_GGT_chr15_73979507  GGT CD276  "$P/GGT_CD276.pdf"
Rscript scripts/figures/figS09_locus_zoom_pg.R locus_GGT_chr16_11644842  GGT LITAF  "$P/GGT_LITAF.pdf"
Rscript scripts/figures/figS09_locus_zoom_pg.R locus_GGT_chr6_53902843   GGT MLIP   "$P/GGT_MLIP.pdf"
Rscript scripts/figures/figS09_locus_zoom_pg.R locus_GGT_chr3_149211512  GGT TM4SF4 "$P/GGT_TM4SF4.pdf"

echo "===== DISEASE / PDFF (EUR) ====="
for g in IL18R1 F13B ERCC2 SHMT1; do
  echo "--- $g ---"
  Rscript scripts/figures/figS09_locus_zoom_disease.R "$g"
done

echo "ALL DONE $(date)"
