#!/usr/bin/env bash
#SBATCH --job-name=masld-enformer-tf217-l40s
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export TF_CPP_MIN_LOG_LEVEL=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/enformer-crested-tf217-l40s-runtime-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Enformer runtime probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

module purge
module load TensorFlow/2.17.0-base
python - <<'PY' > "${STAGE}/runtime_receipt.json"
import json
import os
import platform

import keras
import numpy as np
import tensorflow as tf

os.environ["TF_DETERMINISTIC_OPS"] = "1"
tf.keras.utils.set_random_seed(20260824)
devices = tf.config.list_physical_devices("GPU")
if len(devices) != 1:
    raise SystemExit(f"expected one visible GPU, observed {len(devices)}")
details = tf.config.experimental.get_device_details(devices[0])
name = str(details.get("device_name", ""))
if "L40S" not in name.upper().replace(" ", ""):
    raise SystemExit(f"expected L40S, observed {name!r}")

inputs = tf.reshape(tf.linspace(-1.0, 1.0, 2 * 256 * 4), (2, 256, 4))
layer = keras.layers.Conv1D(16, 15, padding="same", activation="gelu")
with tf.GradientTape() as tape:
    values = layer(inputs, training=True)
    loss = tf.reduce_mean(tf.square(values))
gradients = tape.gradient(loss, layer.trainable_variables)
if any(value is None for value in gradients):
    raise SystemExit("runtime backward produced a missing gradient")
optimizer = keras.optimizers.Adam(1e-4)
optimizer.apply_gradients(zip(gradients, layer.trainable_variables, strict=True))
repeat1 = layer(inputs, training=False).numpy()
repeat2 = layer(inputs, training=False).numpy()
if not np.array_equal(repeat1, repeat2) or not np.isfinite(repeat1).all():
    raise SystemExit("runtime repeat or finiteness check failed")

print(json.dumps({
    "schema_version": "masld-bench-enformer-crested-runtime-probe-v1",
    "status": "pass",
    "python": platform.python_version(),
    "tensorflow": tf.__version__,
    "keras": keras.__version__,
    "gpu_count": len(devices),
    "gpu_name": name,
    "cuda_built": bool(tf.test.is_built_with_cuda()),
    "tiny_forward_shape": list(repeat1.shape),
    "tiny_forward_backward_passed": True,
    "deterministic_repeat_bit_identical": True,
    "model_checkpoint_loaded": False,
    "project_data_read": False,
    "outcomes_read": False,
    "sealed_features_or_labels_read": False,
    "native_sonnet_identity": False,
    "open_champion_eligible": False,
}, sort_keys=True, indent=2))
PY
python --version > "${STAGE}/validation/python.version.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
nvidia-smi -q > "${STAGE}/validation/nvidia-smi.txt"
module list 2> "${STAGE}/validation/module-list.txt"
sha256sum "${ROOT}/slurm/enformer_crested_probe_tf217_l40s.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "runtime_receipt.json").read_text())
if receipt.get("status") != "pass" or receipt.get("model_checkpoint_loaded") is not False:
    raise SystemExit("Enformer CREsted runtime receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "restricted_sequence_model_runtime_probe",
        "model_id": "enformer_crested_restricted_port",
        "tensorflow": receipt["tensorflow"],
        "keras": receipt["keras"],
        "gpu_name": receipt["gpu_name"],
        "model_checkpoint_loaded": False,
        "outcomes_read": False,
        "open_champion_eligible": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Enformer CREsted TF2.17 L40S runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
