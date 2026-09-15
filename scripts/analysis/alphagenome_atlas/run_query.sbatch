#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_query
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/logs/query-%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT=$MASLD_PROJECT_ROOT/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z
export PYTHONNOUSERSITE=1
cd $MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
$PY 06_query_atlas.py direct
echo DIRECT_DONE
$PY 06_query_atlas.py enzyme
echo ENZYME_DONE
