#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-native-probe
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
OVERLAY="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/chrombpnet-v1.0.1-21063022"
SOURCE="${ACQUISITION}/source/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c"
BASE_ARTIFACTS_SHA256=b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229
IMAGE_SHA256=96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa
OVERLAY_ARTIFACTS_SHA256=61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1
SOURCE_ARTIFACTS_SHA256=04bd59c825874f2d2d7dbf38b741633ccacad4735e95997a7477288e9c9656c4
TARGET="${ROOT}/executions/chrombpnet-native-runtime-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PYTHONPATH_VALUE="${OVERLAY}/site-packages:${SOURCE}"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite runtime probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${BASE_ARTIFACTS_SHA256}" "${BASE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${IMAGE_SHA256}" "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' "${OVERLAY_ARTIFACTS_SHA256}" "${OVERLAY}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${ACQUISITION}/ARTIFACTS.json" \
  | sha256sum --check --strict

module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version \
  --format=csv,noheader > "${STAGE}/validation/nvidia-smi.csv"
singularity --version > "${STAGE}/validation/singularity.version.txt"
sha256sum "${ROOT}/slurm/probe_chrombpnet_native_runtime.sbatch" \
  > "${STAGE}/validation/probe-script.sha256"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from __future__ import annotations

from importlib.metadata import version
import hashlib
import json
from pathlib import Path
import platform
from types import SimpleNamespace
import sys

import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp

from chrombpnet.training.models import bpnet_model
from chrombpnet.training.utils.losses import multinomial_nll


output = Path(sys.argv[1]).resolve(strict=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


expected_versions = {
    "absl-py": "1.0.0",
    "gast": "0.5.3",
    "numpy": "1.23.4",
    "protobuf": "3.20.0",
    "six": "1.14.0",
    "tensorflow": "2.8.2",
    "tensorflow-probability": "0.15.0",
}
observed_versions = {name: version(name) for name in expected_versions}
if observed_versions != expected_versions:
    raise RuntimeError(f"runtime versions differ: {observed_versions!r}")

gpus = tf.config.list_physical_devices("GPU")
if len(gpus) != 1:
    raise RuntimeError(f"exactly one GPU is required; observed {gpus!r}")
details = tf.config.experimental.get_device_details(gpus[0])
device_name = str(details.get("device_name", ""))
capability = tuple(details.get("compute_capability", ()))
if "L40S" not in device_name or capability != (8, 9):
    raise RuntimeError(
        f"runtime is not the admitted L40S: {device_name!r}, {capability!r}"
    )

seed = 20260824
tf.keras.backend.clear_session()
tf.random.set_seed(seed)
np.random.seed(seed)
arguments = SimpleNamespace(seed=seed, learning_rate=0.001)
parameters = {
    "filters": 8,
    "n_dil_layers": 2,
    "counts_loss_weight": 1.0,
    "inputlen": 256,
    "outputlen": 128,
}
model = bpnet_model.getModelGivenModelOptionsAndWeightInits(arguments, parameters)
if [tuple(value.shape) for value in model.outputs] != [(None, 128), (None, 1)]:
    raise RuntimeError("native BPNet output shapes differ")

rng = np.random.default_rng(seed)
bases = rng.integers(0, 4, size=(4, 256))
sequence = np.eye(4, dtype=np.float32)[bases]
reverse_complement = sequence[:, ::-1, ::-1].copy()
if not np.array_equal(reverse_complement[:, ::-1, ::-1], sequence):
    raise RuntimeError("reverse-complement orientation fixture is not involutive")
profile_targets = rng.poisson(1.5, size=(4, 128)).astype(np.float32)
profile_targets[:, 0] += 1.0
count_targets = np.log1p(profile_targets.sum(axis=1, keepdims=True)).astype(
    np.float32
)

with tf.device("/GPU:0"):
    first_losses = model.train_on_batch(
        sequence, [profile_targets, count_targets]
    )
    forward = model.predict(sequence, verbose=0)
    reverse_forward = model.predict(reverse_complement, verbose=0)
if not all(np.isfinite(np.asarray(value)).all() for value in first_losses):
    raise RuntimeError("native BPNet forward/backward loss is non-finite")
for values in (*forward, *reverse_forward):
    if not np.isfinite(values).all():
        raise RuntimeError("native BPNet prediction is non-finite")
if [value.shape for value in reverse_forward] != [(4, 128), (4, 1)]:
    raise RuntimeError("reverse-complement forward shapes differ")

checkpoint = output / "checkpoint_after_one_step.h5"
model.save(checkpoint, include_optimizer=True, save_format="h5")
iteration_after_save = int(model.optimizer.iterations.numpy())
restored = tf.keras.models.load_model(
    checkpoint,
    custom_objects={"multinomial_nll": multinomial_nll},
    compile=True,
)
restored_forward = restored.predict(sequence, verbose=0)
reload_max_abs = max(
    float(np.max(np.abs(left - right)))
    for left, right in zip(forward, restored_forward)
)
if reload_max_abs > 1e-6:
    raise RuntimeError(f"checkpoint reload tolerance failed: {reload_max_abs}")
if int(restored.optimizer.iterations.numpy()) != iteration_after_save:
    raise RuntimeError("checkpoint did not restore optimizer iteration")

with tf.device("/GPU:0"):
    original_resume_losses = model.train_on_batch(
        sequence, [profile_targets, count_targets]
    )
    restored_resume_losses = restored.train_on_batch(
        sequence, [profile_targets, count_targets]
    )
    original_resumed = model.predict(sequence, verbose=0)
    restored_resumed = restored.predict(sequence, verbose=0)
resume_max_abs = max(
    float(np.max(np.abs(left - right)))
    for left, right in zip(original_resumed, restored_resumed)
)
if resume_max_abs > 1e-5:
    raise RuntimeError(f"checkpoint resume tolerance failed: {resume_max_abs}")
if int(restored.optimizer.iterations.numpy()) != iteration_after_save + 1:
    raise RuntimeError("restored optimizer did not advance exactly one step")
resumed_checkpoint = output / "checkpoint_after_resume.h5"
restored.save(resumed_checkpoint, include_optimizer=True, save_format="h5")

build_info = tf.sysconfig.get_build_info()
result = {
    "schema_version": "masld-bench-chrombpnet-native-runtime-probe-v1",
    "status": "pass",
    "python": platform.python_version(),
    "versions": observed_versions,
    "tensorflow_probability_import": tfp.__version__,
    "device_name": device_name,
    "device_capability": list(capability),
    "tensorflow_build_info": build_info,
    "source_revision": "eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c",
    "seed": seed,
    "model_parameters": parameters,
    "forward_backward_passed": True,
    "reverse_complement_fixture_passed": True,
    "checkpoint_load_passed": True,
    "checkpoint_resume_passed": True,
    "reload_atol": 1e-6,
    "reload_max_abs": reload_max_abs,
    "resume_atol": 1e-5,
    "resume_max_abs": resume_max_abs,
    "iteration_after_save": iteration_after_save,
    "iteration_after_resume": int(restored.optimizer.iterations.numpy()),
    "first_losses": [float(value) for value in np.atleast_1d(first_losses)],
    "original_resume_losses": [
        float(value) for value in np.atleast_1d(original_resume_losses)
    ],
    "restored_resume_losses": [
        float(value) for value in np.atleast_1d(restored_resume_losses)
    ],
    "checkpoint_sha256": sha256_file(checkpoint),
    "resumed_checkpoint_sha256": sha256_file(resumed_checkpoint),
}
(output / "probe.json").write_text(
    json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
(output / "SOURCE.tsv").write_text(
    "field\tvalue\n"
    "model_id\tchrombpnet\n"
    "release\t1.0.1\n"
    "repository_revision\teaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c\n"
    "base_runtime_artifacts_sha256\t"
    "b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229\n"
    "overlay_artifacts_sha256\t"
    "61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1\n",
    encoding="utf-8",
)
print(json.dumps(result, sort_keys=True))
PY

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
            "artifact_class": "chrombpnet_native_runtime_probe",
            "model_id": "chrombpnet",
            "release": "1.0.1",
            "source_revision": "eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c",
            "device_class": "NVIDIA_L40S_sm89",
            "forward_backward_passed": True,
            "checkpoint_resume_passed": True,
        },
    )
)
PY

mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet native runtime probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
