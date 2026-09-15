#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_pilot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/logs/pilot-%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT=$MASLD_PROJECT_ROOT/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z
export PYTHONNOUSERSITE=1
cd $MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python 05_pilot_atlas.py
echo PILOT_DONE
