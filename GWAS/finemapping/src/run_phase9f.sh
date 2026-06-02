#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=phase9f
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/phase9f_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/phase9f_%j.err
#
# Phase 9f: Aggregate PolyFun SuSiE-COLOC + atlas refresh + 4-way EUR figure regen.

set -eo pipefail
PROJ="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${PROJ}"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

echo "[phase9f] start $(date)"
echo "[phase9f] step 1 — combine PolyFun SuSiE-COLOC CSVs"
Rscript GWAS/finemapping/src/07_combine_susie_coloc_polyfun.R

echo "[phase9f] step 2 — atlas refresh (PolyFun columns)"
Rscript RNA-seq/77b_add_polyfun_atlas_columns.R

echo "[phase9f] step 3 — re-render 4-way EUR LD panel comparison"
Rscript scripts/figures/presentation_4way_eur.R

echo "[phase9f] done $(date)"
