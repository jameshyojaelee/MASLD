#!/usr/bin/env bash
#SBATCH --job-name=masld-sequence-control-unit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONHASHSEED=20260824

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/chrombpnet-v1.0.1-21063022"
SOURCE="${ACQUISITION}/source/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c"
TARGET="${ROOT}/executions/sequence-control-unit-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PYTHONPATH_VALUE="${ROOT}:${CORE}/site-packages:${SOURCE}"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite sequence-control unit artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 04bd59c825874f2d2d7dbf38b741633ccacad4735e95997a7477288e9c9656c4 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict

module load singularity/4.2.2
singularity --version > "${STAGE}/validation/singularity.version.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env PYTHONHASHSEED=20260824 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python -m unittest discover \
    -s "${ROOT}/tests/unit" \
    -p 'sequence_control_test_*.py' \
    -v \
    > "${STAGE}/validation/tests.stdout.txt" \
    2> "${STAGE}/validation/tests.stderr.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/config/sequence_control_architectures.json" \
  "${ROOT}/scripts/sequence_control_architectures.py" \
  "${ROOT}/scripts/sequence_control_prepare_training.py" \
  "${ROOT}/scripts/sequence_control_run_training.py" \
  "${ROOT}/scripts/sequence_control_predict_ccre.py" \
  "${ROOT}/tests/unit/sequence_control_test_models.py" \
  "${ROOT}/slurm/sequence_control_test_models.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
(root / "validation" / "summary.json").write_text(
    json.dumps(
        {
            "schema_version": "masld-bench-sequence-control-unit-v1",
            "status": "pass",
            "models": ["sequence_cnn_control", "sequence_transformer_control"],
            "test_outcomes_used": False,
        },
        sort_keys=True,
        indent=2,
    ) + "\n",
    encoding="utf-8",
)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "sequence_control_unit",
            "models": ["sequence_cnn_control", "sequence_transformer_control"],
            "seed": 20260824,
            "test_outcomes_used": False,
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Sequence-control unit artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
