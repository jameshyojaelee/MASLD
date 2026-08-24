#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-bias-probe
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

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/chrombpnet-v1.0.1-21063022"
SOURCE="${ACQUISITION}/source/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c"
BIAS="${ROOT}/executions/chrombpnet-bias-models-zenodo7443683-21063957"
TARGET="${ROOT}/executions/chrombpnet-pretrained-bias-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PYTHONPATH_VALUE="${CORE}/site-packages:${SOURCE}"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite pretrained-bias probe: %s\n' "${TARGET}" >&2
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
printf '%s  %s\n' 8ecc40f452328df965be823d14f9ce76e2777407dba62cc7826685be635bb6e5 \
  "${BIAS}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 2b69bbbd88ca227ed417d3922cbb2b34c532f4d94da4c01fa8e51e2611520bfc \
  "${ROOT}/executions/chrombpnet-native-runtime-probe-21063590/ARTIFACTS.json" \
  | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BIAS}" \
  > "${STAGE}/validation/bias-tree-verification.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

print(json.dumps(verify_frozen_tree(Path(sys.argv[1])), sort_keys=True))
PY

module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
singularity --version > "${STAGE}/validation/singularity.version.txt"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python - "${BIAS}" "${STAGE}" <<'PY'
from hashlib import sha256
from pathlib import Path
import json
import sys

import numpy as np
import tensorflow as tf

from chrombpnet.training.utils.losses import multinomial_nll


bias_root = Path(sys.argv[1]).resolve(strict=True)
output = Path(sys.argv[2]).resolve(strict=True)
models = {
    "bulk_ATAC_ENCSR868FGK": (
        bias_root / "models/bias_models/ATAC/ENCSR868FGK_bias_fold_0.h5",
        "ae84a033b4c29a502e45e720cd664a0808dc9b98f6d29b35db7158c0a74719aa",
    ),
    "scATAC_dermal_fibroblast": (
        bias_root / "models/bias_models/ATAC/scATAC_dermal_fibroblast.h5",
        "3fcadc521716156405650b1f439bce81c074193c226b4883471b328f1c7b27be",
    ),
}


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


gpus = tf.config.list_physical_devices("GPU")
if len(gpus) != 1:
    raise SystemExit("exactly one GPU is required")
details = tf.config.experimental.get_device_details(gpus[0])
if "L40S" not in str(details.get("device_name")) or tuple(
    details.get("compute_capability", ())
) != (8, 9):
    raise SystemExit(f"unexpected GPU details: {details!r}")

seed = 20260824
results = {}
prediction_vectors = {}
for model_id, (path, expected_sha256) in models.items():
    if sha256_file(path) != expected_sha256:
        raise SystemExit(f"{model_id} checkpoint SHA-256 differs")
    tf.keras.backend.clear_session()
    tf.random.set_seed(seed)
    model = tf.keras.models.load_model(
        path,
        custom_objects={"multinomial_nll": multinomial_nll},
        compile=False,
    )
    input_shape = tuple(model.input_shape)
    if len(input_shape) != 3 or input_shape[0] is not None or input_shape[2] != 4:
        raise SystemExit(f"{model_id} input shape differs: {input_shape!r}")
    input_length = int(input_shape[1])
    model_rng = np.random.default_rng(seed)
    bases = model_rng.integers(0, 4, size=(2, input_length))
    sequence = np.eye(4, dtype=np.float32)[bases]
    reverse_complement = sequence[:, ::-1, ::-1].copy()
    with tf.device("/GPU:0"):
        forward = model.predict(sequence, verbose=0)
        reverse = model.predict(reverse_complement, verbose=0)
    if not isinstance(forward, list) or len(forward) != 2:
        raise SystemExit(f"{model_id} output topology differs")
    output_shapes = [list(value.shape) for value in forward]
    if (
        output_shapes[0][0] != 2
        or output_shapes[1] != [2, 1]
        or any(not np.isfinite(value).all() for value in (*forward, *reverse))
    ):
        raise SystemExit(f"{model_id} output values or shapes differ")
    prediction_vectors[model_id] = np.concatenate(
        [value.reshape(-1) for value in forward]
    )
    results[model_id] = {
        "checkpoint_sha256": expected_sha256,
        "input_shape": list(input_shape),
        "output_shapes": output_shapes,
        "parameters": int(model.count_params()),
        "forward_finite": True,
        "reverse_complement_forward_finite": True,
        "profile_mean": float(np.mean(forward[0])),
        "count_mean": float(np.mean(forward[1])),
    }
left, right = prediction_vectors.values()
if left.shape != right.shape:
    raise SystemExit("ATAC bias checkpoints have incompatible output vectors")
correlation = float(np.corrcoef(left, right)[0, 1])
if not np.isfinite(correlation):
    raise SystemExit("ATAC bias checkpoint correlation is non-finite")
probe = {
    "schema_version": "masld-bench-chrombpnet-pretrained-bias-probe-v1",
    "status": "pass",
    "seed": seed,
    "device_name": str(details.get("device_name")),
    "device_capability": list(details.get("compute_capability", ())),
    "models": results,
    "between_checkpoint_prediction_correlation": correlation,
    "weight_license_status": "unknown_not_redistributable",
    "allowed_role": "restricted_smoke_comparator",
}
(output / "probe.json").write_text(
    json.dumps(probe, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(probe, sort_keys=True))
PY
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum "${ROOT}/slurm/probe_chrombpnet_pretrained_bias_models.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "probe.json").read_text())
print(
    freeze_tree(
        root,
        {
            "artifact_class": "chrombpnet_pretrained_bias_probe",
            "release": "1.0.1",
            "source_revision": "eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c",
            "bias_artifacts_sha256": "8ecc40f452328df965be823d14f9ce76e2777407dba62cc7826685be635bb6e5",
            "models": sorted(probe["models"]),
            "weight_license_status": probe["weight_license_status"],
            "allowed_role": probe["allowed_role"],
            "status": probe["status"],
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet pretrained bias probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
