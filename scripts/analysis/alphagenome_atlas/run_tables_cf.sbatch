#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_tables_CF
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/logs/tablesCF-%j.out
set -euo pipefail
# Packages C (splicing) and F (tissue/ancestry) were in the preview's step list but were left out of
# run_tables.sbatch, so the real run would have shipped without two plan deliverables. They read the same
# archive the earlier steps wrote and add no API calls.
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT=$MASLD_PROJECT_ROOT/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z
export PYTHONNOUSERSITE=1
cd $MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
grep -q TABLES_DONE $AGA_OUT_ROOT/logs/tables-21558781.out
for s in 15_package_C_splicing_join 17_package_F_tissue_ancestry; do
  echo "== $s"; $PY $s.py
done
echo TABLES_CF_DONE
