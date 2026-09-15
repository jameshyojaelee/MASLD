#!/usr/bin/env bash
#SBATCH --job-name=model-test-504
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
runtime=${root}/executions/environments/transcriptformer-torch2.5.1-21070738
runtime_sha256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
printf '%s  %s\n' "${runtime_sha256}" "${runtime}/ARTIFACTS.json" \
  | sha256sum --check --strict
export PYTHONPATH="${root}:${root}/src"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "${root}"
"${runtime}/env/bin/python" -m unittest -v \
  tests.unit.test_fit_predict_common_cell_heads_study_50000
