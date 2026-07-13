#!/bin/bash -l
#SBATCH --job-name=SuSiE
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/runs/uniform35_v2_2026-07-13/logs/aggregate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/runs/uniform35_v2_2026-07-13/logs/aggregate_%j.err

set -euo pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export UNIFORM35_RUN_ROOT="${UNIFORM35_RUN_ROOT:-${PROJECT_ROOT}/GWAS/finemapping/runs/uniform35_v2_2026-07-13}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping
exec Rscript "${PROJECT_ROOT}/GWAS/finemapping/src/seqfunc_v2/finemap/03_aggregate.R"

