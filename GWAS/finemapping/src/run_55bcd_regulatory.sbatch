#!/bin/bash
#SBATCH --job-name=gwasatac
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/run_55bcd_%j.log
# ============================================================================
# run_55bcd_regulatory.sbatch
# Builds the Tier-2 regulatory variant-to-gene CROSS-VALIDATION atlas layers:
#   55b ABC enhancer->gene (downloads Nasser2021 ~340 MB if absent; ~2 GB uncompressed)
#   55c liver caQTL credible-set overlap (Currin 2025)
#   55d SCREEN cCRE overlap (+ HepG2 H3K27ac QC)
#   55e combine -> regulatory_atlas_columns.tsv
# Runs ON A COMPUTE NODE (heavy ABC download/overlap + per-chr caQTL scan).
# ============================================================================
set -euo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC="$BASE/GWAS/finemapping/src"
export MASLD_PROJECT_ROOT="$BASE"

echo "=== node: $(hostname) | start: $(date) ==="
cd "$BASE"

echo ""; echo "######## 55b ABC enhancer->gene ########"
micromamba run -n rnaseq Rscript "$SRC/55b_abc_variant_to_gene.R"

echo ""; echo "######## 55c liver caQTL credible-set ########"
micromamba run -n rnaseq Rscript "$SRC/55c_caqtl_credset.R"

echo ""; echo "######## 55d SCREEN cCRE overlap ########"
micromamba run -n rnaseq Rscript "$SRC/55d_ccre_overlap.R"

echo ""; echo "######## 55e combine -> regulatory_atlas_columns.tsv ########"
micromamba run -n rnaseq Rscript "$SRC/55e_combine_regulatory_atlas.R"

echo ""; echo "=== done: $(date) ==="
