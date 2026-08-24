#!/usr/bin/env bash
#SBATCH --job-name=masld-dedup-tn5-test
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
OVERLAY="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
BASE_ARTIFACTS_SHA256=b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229
IMAGE_SHA256=96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa
OVERLAY_ARTIFACTS_SHA256=61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1

printf '%s  %s\n' "${BASE_ARTIFACTS_SHA256}" "${BASE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${IMAGE_SHA256}" "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' "${OVERLAY_ARTIFACTS_SHA256}" "${OVERLAY}/ARTIFACTS.json" \
  | sha256sum --check --strict
module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${OVERLAY}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m unittest -v tests.unit.test_deduplicated_tn5_bigwigs
