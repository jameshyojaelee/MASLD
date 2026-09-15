#!/usr/bin/env bash
#SBATCH --job-name=model-test-503
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
micromamba=/gpfs/commons/home/jameslee/.local/bin/micromamba
runtime_prefix=/gpfs/commons/home/jameslee/micromamba/envs/immune-atlas
export PYTHONPATH="${root}:${root}/src"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "${root}"
"${micromamba}" run --prefix "${runtime_prefix}" env \
  PYTHONNOUSERSITE=1 PYTHONPATH="${root}:${root}/src" \
  python -m unittest -v tests.unit.test_fit_predict_cell_baselines_study_50000
