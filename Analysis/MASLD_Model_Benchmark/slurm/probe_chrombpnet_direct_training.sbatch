#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-direct-train
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
ACQUISITION="${ROOT}/executions/upstream-acquisitions/chrombpnet-v1.0.1-21063022"
SOURCE="${ACQUISITION}/source/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c"
INPUT="${ROOT}/executions/chrombpnet-train-only-hyperparams-probe-21064020"
TARGET="${ROOT}/executions/chrombpnet-direct-training-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PYTHONPATH_VALUE="${CORE}/site-packages:${SOURCE}"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite direct-training probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/model" "${STAGE}/validation"
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 04bd59c825874f2d2d7dbf38b741633ccacad4735e95997a7477288e9c9656c4 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 5713097b91248632eb8486cd7e6c59d4c1f53832f69f5febca1631d113b480d2 \
  "${INPUT}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${INPUT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
print("input frozen-tree verification passed")
PY

module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
singularity --version > "${STAGE}/validation/singularity.version.txt"
cd "${INPUT}/prepared"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env PYTHONHASHSEED=20260824 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python "${ROOT}/scripts/run_chrombpnet_training.py" \
    --genome "${INPUT}/fixture/fixture.fa" \
    --bigwig "${INPUT}/fixture/fixture.tn5.bw" \
    --peaks train_only_filtered.peaks.bed \
    --nonpeaks train_only_filtered.nonpeaks.bed \
    --output_prefix "${STAGE}/model/chrombpnet" \
    --chr_fold_path actual_fold.json \
    --epochs 2 \
    --early-stop 1 \
    --batch_size 4 \
    --learning-rate 0.001 \
    --params chrombpnet_model_params.actual_fold.tsv \
    --seed 20260824 \
    --architecture_from_file \
      "${SOURCE}/chrombpnet/training/models/chrombpnet_with_bias_model.py" \
    > "${STAGE}/validation/training.stdout.txt" \
    2> "${STAGE}/validation/training.stderr.txt"

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

import numpy as np
import tensorflow as tf

from chrombpnet.training.utils.losses import multinomial_nll

root = Path(sys.argv[1]).resolve(strict=True)
model_path = root / "model" / "chrombpnet.h5"
nobias_path = root / "model" / "chrombpnet_nobias.h5"
if not model_path.is_file() or not nobias_path.is_file():
    raise SystemExit("native training checkpoints are missing")
custom_objects = {"multinomial_nll": multinomial_nll, "tf": tf}
model = tf.keras.models.load_model(model_path, custom_objects=custom_objects, compile=False)
nobias = tf.keras.models.load_model(
    nobias_path, custom_objects=custom_objects, compile=False
)
rng = np.random.default_rng(20260824)
sequence = np.eye(4, dtype=np.float32)[rng.integers(0, 4, size=(2, 2114))]
predictions = model.predict(sequence, verbose=0)
nobias_predictions = nobias.predict(sequence, verbose=0)
if len(predictions) != 2 or len(nobias_predictions) != 2:
    raise SystemExit("native training output topology differs")
if any(not np.isfinite(value).all() for value in predictions + nobias_predictions):
    raise SystemExit("native training emitted non-finite predictions")
with (root / "model" / "chrombpnet.log").open(encoding="utf-8", newline="") as handle:
    history = list(csv.DictReader(handle))
if not 1 <= len(history) <= 2 or any(not np.isfinite(float(row["loss"])) for row in history):
    raise SystemExit("native training history differs")
arguments = json.loads((root / "model" / "chrombpnet.args.json").read_text())
if arguments["epochs"] != 2 or arguments["seed"] != 20260824:
    raise SystemExit("native training argument capture differs")
validation = {
    "schema_version": "masld-bench-chrombpnet-direct-training-probe-v1",
    "status": "pass",
    "epochs_completed": len(history),
    "seed": 20260824,
    "train_outcomes_used": True,
    "validation_outcomes_used_for_early_stopping": True,
    "test_outcomes_used": False,
    "model_forward_passed": True,
    "nobias_model_forward_passed": True,
    "source_bias_role": "restricted_smoke_comparator",
}
(root / "validation" / "probe.json").write_text(
    json.dumps(validation, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(validation, sort_keys=True))
PY
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/run_chrombpnet_training.py" \
  "${ROOT}/slurm/probe_chrombpnet_direct_training.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${INPUT}" \
  "${STAGE}/validation/input-verification.json" <<'PY'
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
            "artifact_class": "chrombpnet_direct_training_probe",
            "model_id": "chrombpnet",
            "release": "1.0.1",
            "source_revision": "eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c",
            "input_artifacts_sha256": "5713097b91248632eb8486cd7e6c59d4c1f53832f69f5febca1631d113b480d2",
            "seed": 20260824,
            "epochs": 2,
            "source_bias_role": "restricted_smoke_comparator",
            "test_outcomes_used": False,
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet direct-training probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
