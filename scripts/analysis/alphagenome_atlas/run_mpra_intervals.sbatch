#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_mpra
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3a-benchmarks-20260909T190214Z/logs/mpra-%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design AGA_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3a-benchmarks-20260909T190214Z PYTHONNOUSERSITE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_atlas
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python 32_mpra_library.py
echo MPRA_DONE
