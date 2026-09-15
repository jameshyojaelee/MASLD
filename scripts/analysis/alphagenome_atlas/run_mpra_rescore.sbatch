#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_mpra_rescore
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3a-benchmarks-20260909T190214Z/logs/mpra_rescore-%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3a-benchmarks-20260909T190214Z
export PYTHONNOUSERSITE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_atlas
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python 33_benchmark_eval.py mpra_rescore
echo MPRA_RESCORE_DONE
