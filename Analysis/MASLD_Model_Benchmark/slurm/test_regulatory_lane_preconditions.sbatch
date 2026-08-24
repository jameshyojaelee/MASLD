#!/usr/bin/env bash
#SBATCH --job-name=masld-reg-preconditions
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
cd "${ROOT}"
printf '%s  %s\n' 94f1fc143d00758cc5b3eae25a31618f1769ca1c183faeb2a2a3dc18ad957576 \
  config/evaluation/bpnet_chrombpnet_native_head_comparison.json \
  | sha256sum --check --strict
printf '%s  %s\n' 6c6486be8157910a2e6afdf07f463309d7602b42e2c6a8a5ce9ac1720df3c38a \
  config/adaptation/corgi.toml | sha256sum --check --strict
printf '%s  %s\n' 796b60caf4acaaf801451a3e91963873d3b6cd22b1892b7bbe367ffe64d263ac \
  config/evaluation/context_backbone_feasibility.json \
  | sha256sum --check --strict
printf '%s  %s\n' c648e18c69762287ea5863c2aef97a9496307def35793d38aa1b9513426db009 \
  tests/contract/test_regulatory_lane_preconditions.py \
  | sha256sum --check --strict

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest \
  tests.contract.test_regulatory_lane_preconditions -v
PYTHONPATH="src" python -m masld_bench.cli registry validate --config-root config
