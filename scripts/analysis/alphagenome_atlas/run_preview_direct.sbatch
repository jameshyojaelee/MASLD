#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_preview
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/preview-direct-only-20260909T191144Z/logs/preview-%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/preview-direct-only-20260909T191144Z
export PYTHONNOUSERSITE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_atlas
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
for s in 07_build_prediction_records 08_posterior_summaries 09_package_A_landscape 10_package_B_chromatin 11_package_D_direction 12_package_G_context 13_catalog_sidecar 15_package_C_splicing_join 17_package_F_tissue_ancestry 14_validate; do
  echo "== $s"; $PY $s.py
done
echo PREVIEW_DONE
