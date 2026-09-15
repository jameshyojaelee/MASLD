#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_saturation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p4-saturation-20260909T191917Z/logs/saturation-%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p4-saturation-20260909T191917Z
export PYTHONNOUSERSITE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_atlas
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python 42_saturation_retrieve.py
echo SATURATION_DONE
