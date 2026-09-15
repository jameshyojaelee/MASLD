#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_tables
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/logs/tables-%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT=$MASLD_PROJECT_ROOT/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z
export PYTHONNOUSERSITE=1
cd $MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
grep -q ENZYME_DONE $AGA_OUT_ROOT/logs/query-21558249.out
$PY 06_query_atlas.py gate
for s in 07_build_prediction_records 08_posterior_summaries 09_package_A_landscape 10_package_B_chromatin 11_package_D_direction 12_package_G_context 13_catalog_sidecar 14_validate; do
  echo "== $s"; $PY $s.py
done
echo TABLES_DONE
