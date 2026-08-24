#!/usr/bin/env bash
#SBATCH --job-name=masld-bpnet-native-probe
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
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
ACQUISITION="${ROOT}/executions/upstream-acquisitions/bpnet-v2.0.0-21064186"
SOURCE="${ACQUISITION}/source/bpnet-f4593eedca51741f25b8c5cc0c0d647faf8c8a3c"
TARGET="${ROOT}/executions/bpnet-native-runtime-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PYTHONPATH_VALUE="${CORE}/site-packages:${SOURCE}"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite BPNet runtime probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/checkpoint" "${STAGE}/validation"

printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 8062574bfaa9067b0aca82ffa3eeae87188140ee5f5842399171ad4535cf88f6 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict

module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
singularity --version > "${STAGE}/validation/singularity.version.txt"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env PYTHONHASHSEED=20260824 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
import random
import sys

import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp

from bpnet.model import arch


root = Path(sys.argv[1]).resolve(strict=True)
seed = 20260824
parameters = {
    "input_len": 2114,
    "output_profile_len": 1000,
    "motif_module_params": {
        "filters": [64],
        "kernel_sizes": [21],
        "padding": "valid",
    },
    "syntax_module_params": {
        "num_dilation_layers": 8,
        "filters": 64,
        "kernel_size": 3,
        "padding": "valid",
        "pre_activation_residual_unit": True,
    },
    "profile_head_params": {"filters": 1, "kernel_size": 75, "padding": "valid"},
    "counts_head_params": {
        "units": [1],
        "dropouts": [0.0],
        "activations": ["linear"],
    },
    "profile_bias_module_params": {"kernel_sizes": [1]},
    "counts_bias_module_params": {},
    "loss_weights": [1.0, 1.0],
    "counts_loss": "MSE",
}
tasks = {
    0: {
        "signal": {"source": ["fixture.tn5.bw"]},
        "loci": {"source": ["fixture.peaks.bed"]},
        "background_loci": {"source": ["fixture.nonpeaks.bed"], "ratio": [1.0]},
        "bias": {"source": [], "smoothing": []},
    }
}


def set_seeds() -> None:
    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)


def build_model():
    set_seeds()
    model = arch.BPNet(tasks, json.loads(json.dumps(parameters)), name_prefix="main")
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.001), loss=None)
    return model


gpus = tf.config.list_physical_devices("GPU")
if len(gpus) != 1:
    raise SystemExit(f"exactly one GPU is required; observed {gpus!r}")
details = tf.config.experimental.get_device_details(gpus[0])
device_name = str(details.get("device_name", ""))
capability = tuple(details.get("compute_capability", ()))
if "L40S" not in device_name or capability != (8, 9):
    raise SystemExit(f"runtime is not L40S sm89: {device_name!r}, {capability!r}")

model = build_model()
if model.input_shape != (None, 2114, 4):
    raise SystemExit(f"BPNet input geometry differs: {model.input_shape!r}")
if [tuple(value.shape) for value in model.outputs] != [
    (None, 1000, 1),
    (None, 1),
]:
    raise SystemExit("BPNet output geometry differs")
if len(model.inputs) != 1 or model.inputs[0].name.split(":")[0] != "sequence":
    raise SystemExit("BPNet unexpectedly accepts control or bias input")
conv_layers = [layer for layer in model.layers if isinstance(layer, tf.keras.layers.Conv1D)]
filters = [layer.filters for layer in conv_layers]
dilations = [layer.dilation_rate[0] for layer in conv_layers]
if filters != [64] * 9 + [1] or dilations != [1, 2, 4, 8, 16, 32, 64, 128, 256, 1]:
    raise SystemExit(f"BPNet convolution census differs: {filters!r}, {dilations!r}")

rng = np.random.default_rng(seed)
sequence = np.eye(4, dtype=np.float32)[rng.integers(0, 4, size=(2, 2114))]
reverse = sequence[:, ::-1, ::-1].copy()
if not np.array_equal(reverse[:, ::-1, ::-1], sequence):
    raise SystemExit("reverse-complement fixture is not involutive")
profile = rng.poisson(1.0, size=(2, 1000, 1)).astype(np.float32)
profile[:, 0, 0] += 1.0
logcounts = np.log1p(profile.sum(axis=1)).astype(np.float32)
targets = {
    "profile_predictions": profile,
    "logcounts_predictions": logcounts,
}
weights = np.ones(2, dtype=np.float32)

first = model.train_on_batch(sequence, targets, sample_weight=weights, return_dict=True)
forward = model.predict(sequence, verbose=0)
reverse_forward = model.predict(reverse, verbose=0)
if any(not np.isfinite(float(value)) for value in first.values()):
    raise SystemExit("BPNet forward/backward loss is non-finite")
if any(not np.isfinite(value).all() for value in (*forward, *reverse_forward)):
    raise SystemExit("BPNet forward output is non-finite")

checkpoint = tf.train.Checkpoint(model=model, optimizer=model.optimizer)
manager = tf.train.CheckpointManager(checkpoint, root / "checkpoint", max_to_keep=2)
checkpoint_path = manager.save(checkpoint_number=1)
iteration_after_save = int(model.optimizer.iterations.numpy())
model.save_weights(root / "checkpoint" / "bpnet.weights.h5")

restored = build_model()
restore_status = tf.train.Checkpoint(
    model=restored, optimizer=restored.optimizer
).restore(checkpoint_path)
restore_status.assert_existing_objects_matched()
restored_forward = restored.predict(sequence, verbose=0)
reload_max_abs = max(
    float(np.max(np.abs(left - right)))
    for left, right in zip(forward, restored_forward)
)
if reload_max_abs > 1e-6:
    raise SystemExit(f"BPNet checkpoint reload tolerance failed: {reload_max_abs}")

model.train_on_batch(sequence, targets, sample_weight=weights, return_dict=True)
restored.train_on_batch(sequence, targets, sample_weight=weights, return_dict=True)
resumed = model.predict(sequence, verbose=0)
restored_resumed = restored.predict(sequence, verbose=0)
resume_max_abs = max(
    float(np.max(np.abs(left - right)))
    for left, right in zip(resumed, restored_resumed)
)
if resume_max_abs > 1e-5:
    raise SystemExit(f"BPNet checkpoint resume tolerance failed: {resume_max_abs}")
if int(restored.optimizer.iterations.numpy()) != iteration_after_save + 1:
    raise SystemExit("BPNet restored optimizer did not advance exactly one step")

cpu_model = build_model()
cpu_model.load_weights(root / "checkpoint" / "bpnet.weights.h5")
with tf.device("/CPU:0"):
    cpu_forward = cpu_model(sequence, training=False)
cross_device_max_abs = max(
    float(np.max(np.abs(left - right.numpy())))
    for left, right in zip(forward, cpu_forward)
)
if cross_device_max_abs > 1e-4:
    raise SystemExit(f"BPNet cross-device tolerance failed: {cross_device_max_abs}")

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


result = {
    "schema_version": "masld-bench-bpnet-native-runtime-probe-v1",
    "status": "pass",
    "source_revision": "f4593eedca51741f25b8c5cc0c0d647faf8c8a3c",
    "declared_upstream_tensorflow": "2.4.1",
    "compatibility_runtime_tensorflow": tf.__version__,
    "compatibility_runtime_tensorflow_probability": tfp.__version__,
    "python": platform.python_version(),
    "device_name": device_name,
    "device_capability": list(capability),
    "seed": seed,
    "input_length": 2114,
    "output_length": 1000,
    "parameter_count": model.count_params(),
    "input_count": len(model.inputs),
    "conv_filters": filters,
    "conv_dilations": dilations,
    "bias_or_control_input": False,
    "forward_backward_passed": True,
    "reverse_complement_fixture_passed": True,
    "checkpoint_load_passed": True,
    "checkpoint_resume_passed": True,
    "reload_max_abs": reload_max_abs,
    "resume_max_abs": resume_max_abs,
    "cross_device_max_abs": cross_device_max_abs,
    "checkpoint_index_sha256": sha256_file(Path(checkpoint_path + ".index")),
    "weights_sha256": sha256_file(root / "checkpoint" / "bpnet.weights.h5"),
}
(root / "validation" / "probe.json").write_text(
    json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
(root / "SOURCE.tsv").write_text(
    "field\tvalue\n"
    "model_id\tbpnet\n"
    "pypi_release\t2.0.0\n"
    "repository_revision\tf4593eedca51741f25b8c5cc0c0d647faf8c8a3c\n"
    "upstream_checkpoint\tnone_project_trained\n"
    "compatibility_runtime\ttensorflow_2.8.2_gpu\n",
    encoding="utf-8",
)
print(json.dumps(result, sort_keys=True))
PY

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum "${ROOT}/slurm/bpnet_probe_native_runtime.sbatch" \
  > "${STAGE}/validation/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(freeze_tree(root, {
    "artifact_class": "bpnet_native_runtime_probe",
    "model_id": "bpnet",
    "release": "2.0.0",
    "source_revision": "f4593eedca51741f25b8c5cc0c0d647faf8c8a3c",
    "compatibility_runtime": "tensorflow_2.8.2_gpu",
    "device_class": "NVIDIA_L40S_sm89",
    "bias_or_control_input": False,
    "forward_backward_passed": True,
    "checkpoint_resume_passed": True,
    "status": "passed",
}))
PY

mv "${STAGE}" "${TARGET}"
printf 'BPNet native runtime probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
