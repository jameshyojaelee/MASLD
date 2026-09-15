#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_packageH
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/logs/packageH-%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_atlas
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python 18_package_H_region_covariates.py   # rnaseq env: the chromatin lanes' env (sklearn via common.py)
echo PACKAGE_H_DONE
