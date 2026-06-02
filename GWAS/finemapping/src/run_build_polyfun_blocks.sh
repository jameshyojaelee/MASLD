#!/bin/bash -l
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=cpu
#SBATCH --array=1-22
#SBATCH --job-name=build_polyfun_blocks
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/build_polyfun_blocks_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/build_polyfun_blocks_%A_%a.err

# Phase 8b — convert PolyFun precomputed UKBB LD (.npz/.gz, 3Mb tiles)
# into per-Berisa-Pickrell-block .bim/.ld files at
# data/ld_ref/polyfun_eur/chr<N>/<bs>.<be>/.

set -eo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${BASE}"

eval "$(micromamba shell hook -s bash)"
# Per CLAUDE.md: Python scripts use `spatial` env (rnaseq has scipy/numpy ABI conflict for Python)
micromamba activate spatial

echo "[run_build_polyfun_blocks] chr${SLURM_ARRAY_TASK_ID}  start $(date)"
python GWAS/finemapping/src/build_polyfun_blocks.py
echo "[run_build_polyfun_blocks] chr${SLURM_ARRAY_TASK_ID}  done  $(date)"
