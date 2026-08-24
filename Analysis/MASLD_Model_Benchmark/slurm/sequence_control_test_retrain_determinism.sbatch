#!/usr/bin/env bash
#SBATCH --job-name=masld-sequence-control-determinism
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
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
FIRST="${ROOT}/executions/sequence-control-training-inference-probe-21064208"
SECOND="${ROOT}/executions/sequence-control-training-inference-probe-21064273"
TARGET="${ROOT}/executions/sequence-control-retrain-determinism-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite determinism artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 9fc50e2fd96285e264c671f603dbcb4e53da44f0966608b0c4f274bc5d958b7b \
  "${FIRST}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' c8b89cbce015d100380f659247a900dba5922bc7d419073e5fd83beb62c822e0 \
  "${SECOND}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${FIRST}" "${SECOND}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("retrain artifacts verified")
PY

module purge
module load singularity/4.2.2
singularity --version > "${STAGE}/validation/singularity.version.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${FIRST}" "${SECOND}" \
  "${STAGE}/validation/summary.json" <<'PY'
from pathlib import Path
import json
import sys

import h5py
import numpy as np

first = Path(sys.argv[1]).resolve(strict=True)
second = Path(sys.argv[2]).resolve(strict=True)


def h5_datasets(path: Path):
    values = {}
    with h5py.File(path, "r") as handle:
        handle.visititems(
            lambda name, item: values.__setitem__(name, item[()])
            if isinstance(item, h5py.Dataset)
            else None
        )
    return values


results = {}
for model_id in ("sequence_cnn_control", "sequence_transformer_control"):
    first_model = h5_datasets(first / "models" / f"{model_id}.h5")
    second_model = h5_datasets(second / "models" / f"{model_id}.h5")
    if set(first_model) != set(second_model):
        raise SystemExit(f"{model_id} trained HDF5 dataset roster differs")
    maximum_weight_difference = 0.0
    for key in first_model:
        left, right = first_model[key], second_model[key]
        if left.shape != right.shape or left.dtype != right.dtype:
            raise SystemExit(f"{model_id} trained HDF5 tensor schema differs")
        if np.issubdtype(left.dtype, np.number):
            if left.size:
                maximum_weight_difference = max(
                    maximum_weight_difference,
                    float(np.max(np.abs(left.astype(float) - right.astype(float)))),
                )
            if not np.array_equal(left, right):
                raise SystemExit(f"{model_id} trained tensors are not bit-identical")
        elif not np.array_equal(left, right):
            raise SystemExit(f"{model_id} trained metadata tensors differ")
    for relative in (
        Path("models") / f"{model_id}.log",
        Path("predictions") / model_id / "regional_counts.tsv",
        Path("predictions") / model_id / "standardized_predictions.tsv",
        Path("predictions") / model_id / "summary.json",
        Path("predictions") / model_id / "prediction_contract.json",
    ):
        if (first / relative).read_bytes() != (second / relative).read_bytes():
            raise SystemExit(f"{model_id} retrain artifact differs: {relative}")
    first_profiles = h5_datasets(
        first / "predictions" / model_id / "profile_probabilities.h5"
    )
    second_profiles = h5_datasets(
        second / "predictions" / model_id / "profile_probabilities.h5"
    )
    if set(first_profiles) != set(second_profiles) or any(
        not np.array_equal(first_profiles[key], second_profiles[key])
        for key in first_profiles
    ):
        raise SystemExit(f"{model_id} retrain profiles are not bit-identical")
    results[model_id] = {
        "trained_hdf5_datasets": len(first_model),
        "maximum_weight_absolute_difference": maximum_weight_difference,
        "regional_counts_bit_identical": True,
        "standardized_predictions_bit_identical": True,
        "profile_datasets_bit_identical": True,
    }
summary = {
    "schema_version": "masld-bench-sequence-control-retrain-determinism-v1",
    "status": "pass",
    "seed": 20260824,
    "same_l40s_runtime": True,
    "models": results,
    "trained_predictions_bit_identical": True,
    "held_donor_atac_exposed": False,
    "test_outcomes_used": False,
}
Path(sys.argv[3]).write_text(
    json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(summary, sort_keys=True))
PY
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum "${ROOT}/slurm/sequence_control_test_retrain_determinism.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "sequence_control_retrain_determinism",
            "models": ["sequence_cnn_control", "sequence_transformer_control"],
            "seed": 20260824,
            "trained_predictions_bit_identical": True,
            "test_outcomes_used": False,
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Sequence-control retrain determinism: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
