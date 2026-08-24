#!/usr/bin/env bash
#SBATCH --job-name=masld-rna-atac-reg-eval-runtime
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
EVALUATOR="${ROOT}/src/masld_bench/evaluators/rna_atac_regularized_development.py"
EVALUATOR_SHA256=d824670898285794ddeac92ec19b08b40b56bc37f90fc1c0ce6b728ff1874981

cd "${ROOT}"
printf '%s  %s\n' "${EVALUATOR_SHA256}" "${EVALUATOR}" | sha256sum -c -
/gpfs/commons/home/jameslee/.local/bin/micromamba run \
  --prefix /gpfs/commons/home/jameslee/micromamba/envs/immune-atlas \
  env \
  PYTHONNOUSERSITE=1 \
  python "${EVALUATOR}" --help
