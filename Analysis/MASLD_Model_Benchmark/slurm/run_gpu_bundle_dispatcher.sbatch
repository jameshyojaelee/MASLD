#!/usr/bin/env bash
#SBATCH --job-name=model-dispatcher
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
QUEUE="${ROOT}/config/campaigns/gpu_bundle_queue"
STATE="${ROOT}/executions/gpu-bundle-dispatch-state"
mkdir -p "${QUEUE}" "${STATE}"
module load python/3.11.5-GCCcore-13.2.0
POLICY="${ROOT}/config/campaigns/gpu_dependency_dispatcher_policy_20260825.json"
PYTHONPATH="src:." python -m unittest \
  tests.unit.test_gpu_dependency_dispatcher \
  tests.unit.test_gpu_bundle_dispatcher
sha256sum \
  "${ROOT}/scripts/gpu_dependency_dispatcher.py" \
  "${ROOT}/scripts/gpu_bundle_dispatcher.py" \
  "${ROOT}/tests/unit/test_gpu_dependency_dispatcher.py" \
  "${ROOT}/tests/unit/test_gpu_bundle_dispatcher.py" \
  "${ROOT}/slurm/run_gpu_bundle_dispatcher.sbatch" \
  "${POLICY}"
PYTHONPATH="src:." python "${ROOT}/scripts/gpu_dependency_dispatcher.py" \
  --root "${ROOT}" --queue "${QUEUE}" --state "${STATE}" \
  --policy "${POLICY}" --user "${USER}"
