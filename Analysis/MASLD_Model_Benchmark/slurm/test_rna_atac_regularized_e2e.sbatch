#!/usr/bin/env bash
#SBATCH --job-name=masld-rna-atac-reg-e2e
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark

/gpfs/commons/home/jameslee/.local/bin/micromamba run \
  --prefix /gpfs/commons/home/jameslee/micromamba/envs/immune-atlas \
  env \
  PYTHONNOUSERSITE=1 \
  OMP_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 \
  python tests/scientific/test_rna_atac_regularized_e2e.py \
    --environment-lock executions/environment-locks/cpu_baseline_smoke--8a7dfd2a15a635d5/environment.lock.json \
    --data executions/gse296875-rna-atac-smoke-20260824T015422Z-21039865/subset/gse296875_rna_atac_smoke_1000.h5 \
    --output executions/tests/rna-atac-regularized-e2e-${SLURM_JOB_ID}
