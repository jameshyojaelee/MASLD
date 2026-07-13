#!/usr/bin/env bash
#SBATCH --job-name=MPRA
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --array=0-1
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/MPRA_source_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/MPRA_source_%A_%a.err
set -euo pipefail
ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
cells=(HepG2 LX2)
export MASLD_PROJECT_ROOT=${ROOT}
export MPRA_CELL_LINE=${cells[$SLURM_ARRAY_TASK_ID]}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
micromamba run -n rnaseq Rscript "${ROOT}/GWAS/finemapping/src/seqfunc_v2/mpra/02_source_reproduction.R"
