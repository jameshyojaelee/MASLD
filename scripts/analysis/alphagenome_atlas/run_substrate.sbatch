#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_substrate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/logs/substrate-%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT=$MASLD_PROJECT_ROOT/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z
export PYTHONNOUSERSITE=1
cd $MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
$PY 01_build_universes.py
$PY 02_extract_variants.py
$RS 03_liftover_universeB.R
$PY 04_build_variant_crosswalk.py
echo SUBSTRATE_DONE
