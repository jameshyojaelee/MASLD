#!/usr/bin/env bash
#SBATCH --job-name=masld-bpnet-contract-test
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"

printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict

module load singularity/4.2.2
cd "${ROOT}"
bash -n \
  slurm/bpnet_acquire_v2_0_0.sbatch \
  slurm/bpnet_probe_native_runtime.sbatch \
  slurm/bpnet_probe_training_inference.sbatch \
  slurm/bpnet_probe_training_resume.sbatch \
  slurm/bpnet_test_contract.sbatch \
  slurm/bpnet_train_predict_hepatocyte_donor0_genomic0.sbatch
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${CORE}/site-packages:${ROOT}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m unittest -v \
    tests.unit.test_bpnet_contract \
    tests.unit.test_bpnet_predict_ccre_windows

singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${CORE}/site-packages:${ROOT}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m py_compile \
    scripts/bpnet_contract.py \
    scripts/bpnet_run_training.py \
    scripts/bpnet_predict_ccre_windows.py
