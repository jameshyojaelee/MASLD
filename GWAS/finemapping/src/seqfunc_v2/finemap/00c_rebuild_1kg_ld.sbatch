#!/bin/bash
#SBATCH --job-name=plink2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/ld_rebuild/rebuild_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/ld_rebuild/rebuild_%A_%a.err

# One 1000G LD block per array task.  Output goes to a PARALLEL panel root
# `1kg_<pop>_gram/`; the originals under `1kg_<pop>/` are shared with the
# canonical COLOC pipeline and are never touched.
#
# 128G covers the largest block (M = 21,439 -> a 3.7 GB dense matrix plus
# intermediates).  8 CPUs because cor() is BLAS-3 and threads well.
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
export MASLD_PROJECT_ROOT
export UNIFORM35_RUN_ROOT="${UNIFORM35_RUN_ROOT:-${PROJECT_ROOT}/GWAS/finemapping/runs/uniform35_v3_2026-08-03}"
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "[$(date -u +%FT%TZ)] rebuilding block row ${SLURM_ARRAY_TASK_ID}"
Rscript "${PROJECT_ROOT}/GWAS/finemapping/src/seqfunc_v2/finemap/00c_rebuild_1kg_ld.R" "${SLURM_ARRAY_TASK_ID}"
echo "[$(date -u +%FT%TZ)] done"
