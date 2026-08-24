#!/usr/bin/env bash
#SBATCH --job-name=masld-sequence-control-b32
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
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
SMOKE="${ROOT}/executions/sequence-control-training-inference-probe-21064208"
TARGET="${ROOT}/executions/sequence-control-batch32-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PYTHONPATH_VALUE="${ROOT}:${CORE}/site-packages:${SOURCE}"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite batch-32 probe: %s\n' "${TARGET}" >&2
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
printf '%s  %s\n' 9fc50e2fd96285e264c671f603dbcb4e53da44f0966608b0c4f274bc5d958b7b \
  "${SMOKE}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SMOKE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
print("smoke frozen-tree verification passed")
PY

module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total \
  --format=csv,noheader > "${STAGE}/validation/nvidia-smi.before.csv"
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
  "${IMAGE}" python - "${SMOKE}" "${STAGE}/validation/summary.json" <<'PY'
from pathlib import Path
import json
import random
import subprocess
import sys

import numpy as np
import tensorflow as tf

from chrombpnet.training.utils.losses import multinomial_nll

smoke = Path(sys.argv[1]).resolve(strict=True)
rng = np.random.default_rng(20260824)
sequence = np.eye(4, dtype=np.float32)[rng.integers(0, 4, size=(32, 2114))]
profile = rng.poisson(2.0, size=(32, 1000)).astype(np.float32)
logcount = np.log1p(profile.sum(axis=1, keepdims=True)).astype(np.float32)
results = {}
for model_id in ("sequence_cnn_control", "sequence_transformer_control"):
    random.seed(20260824)
    np.random.seed(20260824)
    tf.random.set_seed(20260824)
    model = tf.keras.models.load_model(
        smoke / "models" / f"{model_id}.h5",
        custom_objects={"multinomial_nll": multinomial_nll, "tf": tf},
        compile=False,
    )
    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
        loss=[multinomial_nll, "mse"],
        loss_weights=[1.0, 210.0],
    )
    loss = np.asarray(model.train_on_batch(sequence, [profile, logcount]))
    if not np.isfinite(loss).all():
        raise SystemExit(f"{model_id} batch-32 backward emitted non-finite loss")
    memory_used = subprocess.check_output(
        [
            "nvidia-smi",
            "--query-gpu=memory.used",
            "--format=csv,noheader,nounits",
        ],
        text=True,
    ).strip()
    results[model_id] = {
        "loss": [float(value) for value in loss.reshape(-1)],
        "gpu_memory_used_mib_after_step": int(memory_used),
    }
    del model
    tf.keras.backend.clear_session()
summary = {
    "schema_version": "masld-bench-sequence-control-batch32-probe-v1",
    "status": "pass",
    "batch_size": 32,
    "input_length": 2114,
    "output_length": 1000,
    "forward_backward_passed": True,
    "synthetic_outcomes_only": True,
    "held_donor_atac_exposed": False,
    "test_outcomes_used": False,
    "models": results,
}
Path(sys.argv[2]).write_text(
    json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(summary, sort_keys=True))
PY
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum "${ROOT}/slurm/sequence_control_probe_batch32.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${SMOKE}" "${STAGE}/validation/input-verification.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

Path(sys.argv[2]).write_text(
    json.dumps(verify_frozen_tree(Path(sys.argv[1])), sort_keys=True, indent=2) + "\n",
    encoding="utf-8",
)
PY
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "sequence_control_batch32_probe",
            "models": ["sequence_cnn_control", "sequence_transformer_control"],
            "batch_size": 32,
            "seed": 20260824,
            "synthetic_outcomes_only": True,
            "test_outcomes_used": False,
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Sequence-control batch-32 probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
