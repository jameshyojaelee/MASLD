#!/usr/bin/env bash
#SBATCH --job-name=squidpy
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/hmsma_class_spatial_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/hmsma_class_spatial_%j.out
set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Reuses the stored per-gene Moran's I values from the 2026-08-10 run, so no
# squidpy recomputation is needed. Only the inference stage changed: the invalid
# across-array t-test is replaced by a matched-gene-set null. Drop
# --from-per-gene to recompute Moran's I from the h5ad objects instead.
micromamba run -n spatial python \
  Analysis/Spatial/scripts/55_hmsma_class_spatial_autocorr.py \
  --from-per-gene Analysis/Spatial/results/hmsma_class_spatial/per_gene_morans_i.tsv.gz \
  --n-null 9999 \
  --seed 42
echo "EXIT=$?"
