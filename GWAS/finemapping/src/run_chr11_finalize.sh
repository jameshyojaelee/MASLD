#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=pf11fin
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/pf11fin_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/pf11fin_%j.err
#
# After 16 chr11 PolyFun workers finish:
#   1. Merge worker_chr11_w*.csv → susie_coloc_chr11.csv
#   2. Re-run Phase 9f aggregator (07_combine + 77b atlas + 4-way EUR figure)

set -eo pipefail
PROJ="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${PROJ}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

echo "[pf11fin] start $(date)"
echo "[pf11fin] step 1 — merge chr11 worker outputs"
Rscript GWAS/finemapping/src/merge_polyfun_chr11_workers.R

echo "[pf11fin] step 2 — re-run Phase 9f (combine + atlas + 4-way EUR figure)"
Rscript GWAS/finemapping/src/07_combine_susie_coloc_polyfun.R
Rscript RNA-seq/77b_add_polyfun_atlas_columns.R
Rscript scripts/figures/presentation_4way_eur.R

echo "[pf11fin] done $(date)"
