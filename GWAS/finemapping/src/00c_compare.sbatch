#!/bin/bash
# Overlap/difference comparison between the Sveinbjornsson-2022 bundle and our existing GWAS.
#SBATCH --job-name=gwas_overlap
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/gwas_overlap_%j.out
set -uo pipefail
SRC=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/gwas-review/GWAS/finemapping/src
micromamba run -n rnaseq Rscript "$SRC/00c_compare_overlapping_gwas.R"
