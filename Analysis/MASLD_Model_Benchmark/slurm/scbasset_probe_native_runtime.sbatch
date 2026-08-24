#!/usr/bin/env bash
#SBATCH --job-name=masld-scbasset-native-probe
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
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/scbasset-aed3a6f713091fd988196b297e9c06e092ff1d22-21064184"
SOURCE="${ACQUISITION}/source/scBasset-aed3a6f713091fd988196b297e9c06e092ff1d22"
TARGET="${ROOT}/executions/scbasset-native-runtime-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
BASE_ARTIFACTS_SHA256=b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229
IMAGE_SHA256=96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa
CORE_ARTIFACTS_SHA256=61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1
SOURCE_ARTIFACTS_SHA256=9b0b418cea67bc788d7657913965c5e11ed62e36c35a8104c381c9007f6dffc5

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite scBasset runtime probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/checkpoints"

printf '%s  %s\n' "${BASE_ARTIFACTS_SHA256}" "${BASE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${IMAGE_SHA256}" "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' "${CORE_ARTIFACTS_SHA256}" "${CORE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${ACQUISITION}/ARTIFACTS.json" \
  | sha256sum --check --strict

module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
singularity --version > "${STAGE}/validation/singularity.version.txt"
sha256sum "${ROOT}/slurm/scbasset_probe_native_runtime.sbatch" \
  > "${STAGE}/validation/probe-script.sha256"

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${CORE}/site-packages" \
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
  --env "PYTHONPATH=${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python - "${SOURCE}" "${STAGE}" <<'PY'
from __future__ import annotations

import ast
from hashlib import sha256
import importlib.util
from importlib.metadata import version
import json
from pathlib import Path
import platform
import sys
from types import ModuleType

import h5py
import numpy as np
import tensorflow as tf


source = Path(sys.argv[1]).resolve(strict=True)
output = Path(sys.argv[2]).resolve(strict=True)


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


expected_versions = {
    "h5py": "3.6.0",
    "numpy": "1.23.4",
    "tensorflow": "2.8.2",
}
observed_versions = {name: version(name) for name in expected_versions}
if observed_versions != expected_versions:
    raise RuntimeError(f"scBasset runtime versions differ: {observed_versions!r}")

gpus = tf.config.list_physical_devices("GPU")
if len(gpus) != 1:
    raise RuntimeError(f"exactly one GPU is required: {gpus!r}")
details = tf.config.experimental.get_device_details(gpus[0])
device_name = str(details.get("device_name", ""))
capability = tuple(details.get("compute_capability", ()))
if "L40S" not in device_name or capability != (8, 9):
    raise RuntimeError(
        f"runtime is not the admitted L40S: {device_name!r}, {capability!r}"
    )

utils_path = source / "scbasset" / "utils.py"
basenji_path = source / "scbasset" / "basenji_utils.py"
if sha256_file(utils_path) != "7b653646fc83867cf38c51ae14fdcdf523f0dc99685d1b9446ea5940de71543c":
    raise RuntimeError("registered scBasset architecture source changed")
if sha256_file(basenji_path) != "165fdc3c2dd85afd300bf71e27199ae907b005e325c33fa0cc1e9b96be9a92b2":
    raise RuntimeError("registered Basenji helper source changed")
tree = ast.parse(utils_path.read_text(encoding="utf-8"))
make_model_nodes = [
    node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "make_model"
]
if len(make_model_nodes) != 1:
    raise RuntimeError("registered scBasset make_model function is not unique")
node = make_model_nodes[0]
if [argument.arg for argument in node.args.args] != [
    "bottleneck_size", "n_cells", "seq_len", "show_summary"
]:
    raise RuntimeError("registered scBasset make_model signature changed")

# The architecture helper imports pysam only for optional FASTA extraction.
# The native runtime probe exercises no FASTA I/O, so a nonfunctional module
# placeholder avoids adding a production dependency while preserving the exact
# audited helper implementation.
sys.modules["pysam"] = ModuleType("pysam")
spec = importlib.util.spec_from_file_location("scbasset_basenji_pinned", basenji_path)
if spec is None or spec.loader is None:
    raise RuntimeError("cannot load pinned scBasset Basenji helper")
basenji = importlib.util.module_from_spec(spec)
spec.loader.exec_module(basenji)


def make_native_model(bottleneck_size: int, n_cells: int) -> tf.keras.Model:
    sequence = tf.keras.Input(shape=(1344, 4), name="sequence")
    current, reverse_bool = basenji.StochasticReverseComplement()(sequence)
    current = basenji.StochasticShift(3)(current)
    current = basenji.conv_block(
        current, filters=288, kernel_size=17, pool_size=3
    )
    current = basenji.conv_tower(
        current,
        filters_init=288,
        filters_mult=1.122,
        repeat=6,
        kernel_size=5,
        pool_size=2,
    )
    current = basenji.conv_block(current, filters=256, kernel_size=1)
    current = basenji.dense_block(
        current, flatten=True, units=bottleneck_size, dropout=0.2
    )
    current = basenji.GELU()(current)
    current = basenji.final(current, units=n_cells, activation="sigmoid")
    current = basenji.SwitchReverse()([current, reverse_bool])
    current = tf.keras.layers.Flatten()(current)
    return tf.keras.Model(inputs=sequence, outputs=current)


seed = 20260824
np.random.seed(seed)
tf.random.set_seed(seed)
tf.keras.backend.clear_session()
model = make_native_model(32, 10)
if tuple(model.input_shape) != (None, 1344, 4) or tuple(model.output_shape) != (None, 10):
    raise RuntimeError("native scBasset input/output shape changed")
model.compile(
    loss=tf.keras.losses.BinaryCrossentropy(),
    optimizer=tf.keras.optimizers.Adam(
        learning_rate=0.01, beta_1=0.95, beta_2=0.9995
    ),
)

rng = np.random.default_rng(seed)
bases = rng.integers(0, 4, size=(4, 1344))
sequences = np.eye(4, dtype=np.float32)[bases]
reverse_complements = sequences[:, ::-1, ::-1].copy()
targets = rng.integers(0, 2, size=(4, 10)).astype(np.float32)
targets[:, 0] = np.array([0, 1, 0, 1], dtype=np.float32)
with tf.device("/GPU:0"):
    first_loss = float(model.train_on_batch(sequences, targets))
    forward = model.predict(sequences, batch_size=2, verbose=0)
    reverse_forward = model.predict(reverse_complements, batch_size=2, verbose=0)
if not np.isfinite(first_loss) or not np.isfinite(forward).all() or not np.isfinite(reverse_forward).all():
    raise RuntimeError("native scBasset forward/backward is non-finite")
if forward.shape != (4, 10) or np.any(forward <= 0) or np.any(forward >= 1):
    raise RuntimeError("native scBasset sigmoid prediction contract changed")
rc_tta = 0.5 * (forward + reverse_forward)
if not np.isfinite(rc_tta).all():
    raise RuntimeError("reverse-complement test-time average is non-finite")

weights_path = output / "checkpoints" / "native_after_one_step.weights.h5"
model.save_weights(weights_path)
dense_candidates = [
    layer
    for layer in model.layers
    if isinstance(layer, tf.keras.layers.Dense)
    and layer.units == 10
    and layer.activation == tf.keras.activations.sigmoid
]
if len(dense_candidates) != 1:
    raise RuntimeError("training-cell output dense layer is not unique")
kernel, intercept = dense_candidates[0].get_weights()
if kernel.shape != (32, 10) or intercept.shape != (10,):
    raise RuntimeError("transductive cell output weight shape changed")
cell_embedding = kernel.T
if cell_embedding.shape != (10, 32):
    raise RuntimeError("training-cell embedding shape changed")

h5_members = []
with h5py.File(weights_path, "r") as handle:
    def record(name: str, value: object) -> None:
        if isinstance(value, h5py.Dataset):
            array = np.asarray(value[()])
            h5_members.append(
                {
                    "name": name,
                    "shape": list(array.shape),
                    "dtype": str(array.dtype),
                    "sha256": sha256(array.tobytes(order="C")).hexdigest(),
                }
            )
    handle.visititems(record)
if not h5_members:
    raise RuntimeError("local HDF5 checkpoint has no tensor datasets")
(output / "checkpoints" / "native_after_one_step.h5_members.json").write_text(
    json.dumps(h5_members, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)

tf.keras.backend.clear_session()
np.random.seed(seed + 1)
tf.random.set_seed(seed + 1)
restored = make_native_model(32, 10)
restored.compile(
    loss=tf.keras.losses.BinaryCrossentropy(),
    optimizer=tf.keras.optimizers.Adam(
        learning_rate=0.01, beta_1=0.95, beta_2=0.9995
    ),
)
restored.load_weights(weights_path)
restored_forward = restored.predict(sequences, batch_size=2, verbose=0)
reload_max_abs = float(np.max(np.abs(forward - restored_forward)))
if reload_max_abs > 1e-6:
    raise RuntimeError(f"scBasset HDF5 reload tolerance failed: {reload_max_abs}")
with tf.device("/GPU:0"):
    resume_loss = float(restored.train_on_batch(sequences, targets))
if not np.isfinite(resume_loss) or int(restored.optimizer.iterations.numpy()) != 1:
    raise RuntimeError("scBasset checkpoint resume step failed")
resumed_path = output / "checkpoints" / "native_after_resume.weights.h5"
restored.save_weights(resumed_path)

result = {
    "schema_version": "masld-bench-scbasset-native-runtime-probe-v1",
    "status": "pass",
    "model_id": "scbasset",
    "source_revision": "aed3a6f713091fd988196b297e9c06e092ff1d22",
    "source_utils_sha256": sha256_file(utils_path),
    "source_basenji_helpers_sha256": sha256_file(basenji_path),
    "python": platform.python_version(),
    "versions": observed_versions,
    "device_name": device_name,
    "device_capability": list(capability),
    "seed": seed,
    "sequence_length": 1344,
    "bottleneck_dimension": 32,
    "training_cell_outputs": 10,
    "parameter_count": int(model.count_params()),
    "first_loss": first_loss,
    "resume_loss": resume_loss,
    "forward_backward_passed": True,
    "reverse_complement_tta_passed": True,
    "hdf5_weights_only_load_passed": True,
    "checkpoint_resume_step_passed": True,
    "reload_atol": 1e-6,
    "reload_max_abs": reload_max_abs,
    "hdf5_tensor_members": len(h5_members),
    "weights_sha256": sha256_file(weights_path),
    "resumed_weights_sha256": sha256_file(resumed_path),
    "cell_embedding_interpretation": "transposed_final_dense_kernel_for_training_cells_only",
    "held_cell_embedding_available": False,
    "ambiguous_sequence_policy": "reject_non_ACGT_before_one_hot_encoding",
    "public_tutorial_weights_used": False,
    "selection_objective_for_campaign": "minimum_binary_crossentropy_on_genomic_valid_ccres_using_nested_training_donor_cell_outputs",
}
(output / "probe.json").write_text(
    json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
(output / "SOURCE.tsv").write_text(
    "field\tvalue\n"
    "model_id\tscbasset\n"
    "repository\tcalico/scBasset\n"
    "repository_revision\taed3a6f713091fd988196b297e9c06e092ff1d22\n"
    "code_license\tApache-2.0\n"
    "training_initialization\tfrom_scratch\n"
    "public_tutorial_weights_used\tfalse\n",
    encoding="utf-8",
)
print(json.dumps(result, sort_keys=True))
PY

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
            "artifact_class": "scbasset_native_runtime_probe",
            "model_id": "scbasset",
            "source_revision": "aed3a6f713091fd988196b297e9c06e092ff1d22",
            "device_class": "NVIDIA_L40S_sm89",
            "forward_backward_passed": True,
            "hdf5_weights_only_load_passed": True,
            "checkpoint_resume_step_passed": True,
            "public_tutorial_weights_used": False,
        },
    )
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scBasset native runtime probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
