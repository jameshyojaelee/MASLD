#!/usr/bin/env bash
#SBATCH --job-name=MPRAnalyze
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --time=90:00:00
#SBATCH --array=0-1
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/MPRAnalyze_fit_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/MPRAnalyze_fit_%A_%a.err
set -euo pipefail
ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
cells=(HepG2 LX2)
export MASLD_PROJECT_ROOT=${ROOT}
export MPRA_CELL_LINE=${cells[$SLURM_ARRAY_TASK_ID]}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
micromamba run -p "${ROOT}/.mamba/seqfunc_mpra_v2" Rscript \
  "${ROOT}/GWAS/finemapping/src/seqfunc_v2/mpra/03_mpranalyze.R"
