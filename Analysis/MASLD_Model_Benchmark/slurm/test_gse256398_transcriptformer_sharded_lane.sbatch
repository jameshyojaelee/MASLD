#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-test-257
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
cd "${ROOT}"
export PYTHONPATH=${ROOT}:${ROOT}/src
"${RUNTIME}/env/bin/python" -m unittest -v \
  tests.unit.test_lock_gse256398_transcriptformer_membership \
  tests.unit.test_build_gse256398_transcriptformer_shards \
  tests.unit.test_transcriptformer_extract_gse256398_shards
"${RUNTIME}/env/bin/python" -m py_compile \
  scripts/lock_gse256398_transcriptformer_membership.py \
  scripts/build_gse256398_transcriptformer_shards.py \
  scripts/transcriptformer_extract_gse256398_shards.py
git diff --check -- \
  scripts/lock_gse256398_transcriptformer_membership.py \
  scripts/build_gse256398_transcriptformer_shards.py \
  scripts/transcriptformer_extract_gse256398_shards.py \
  tests/unit/test_lock_gse256398_transcriptformer_membership.py \
  tests/unit/test_build_gse256398_transcriptformer_shards.py \
  tests/unit/test_transcriptformer_extract_gse256398_shards.py \
  slurm/lock_gse256398_transcriptformer_membership.sbatch \
  slurm/build_gse256398_transcriptformer_shards.sbatch \
  slurm/run_transcriptformer_gse256398_shards.sbatch \
  slurm/test_gse256398_transcriptformer_sharded_lane.sbatch
