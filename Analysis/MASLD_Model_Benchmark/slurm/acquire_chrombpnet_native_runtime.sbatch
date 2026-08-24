#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-runtime
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
IMAGE_URI=docker://tensorflow/tensorflow:2.8.2-gpu
TARGET="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
IMAGE="${STAGE}/tensorflow-tensorflow-2.8.2-gpu.sif"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite runtime output: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/runtime" "${STAGE}/validation"

module load singularity/4.2.2
module load python/3.11.5-GCCcore-13.2.0
export SINGULARITY_CACHEDIR="${TMPDIR:-/tmp}/singularity-cache-${SLURM_JOB_ID}"
export SINGULARITY_TMPDIR="${TMPDIR:-/tmp}"

singularity pull --disable-cache "${IMAGE}" "${IMAGE_URI}"
test -s "${IMAGE}"
sha256sum "${IMAGE}" > "${STAGE}/validation/image.sha256"
singularity --version > "${STAGE}/validation/singularity.version.txt"
singularity inspect --json "${IMAGE}" > "${STAGE}/validation/image.inspect.json"
singularity exec --cleanenv "${IMAGE}" python - <<'PY' \
  > "${STAGE}/validation/runtime_probe.json"
import json
import platform
import sys

import google.protobuf
import h5py
import numpy
import tensorflow as tf

if tf.__version__ != "2.8.2":
    raise SystemExit(f"unexpected TensorFlow version: {tf.__version__}")
result = {
    "gpu_visible_during_cpu_acquisition": bool(tf.config.list_physical_devices("GPU")),
    "h5py_version": h5py.__version__,
    "numpy_version": numpy.__version__,
    "platform": platform.platform(),
    "protobuf_version": google.protobuf.__version__,
    "python_version": sys.version,
    "tensorflow_git_version": getattr(tf.version, "GIT_VERSION", "unknown"),
    "tensorflow_version": tf.__version__,
}
print(json.dumps(result, sort_keys=True, separators=(",", ":")))
PY
singularity exec --cleanenv "${IMAGE}" python -m pip freeze \
  > "${STAGE}/runtime/pip-freeze.txt"

cat > "${STAGE}/SOURCE.tsv" <<EOF
field\tvalue
runtime_id\tchrombpnet_native_tf2_8_2
image_uri\t${IMAGE_URI}
upstream_project\ttensorflow/tensorflow
upstream_version\t2.8.2-gpu
upstream_code_license\tApache-2.0
intended_model\tchrombpnet_1.0.1
installed_chrombpnet\tfalse
EOF

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation" / "runtime_probe.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "chrombpnet_native_runtime_acquisition",
        "image_sha256": sha256_file(
            root / "tensorflow-tensorflow-2.8.2-gpu.sif"
        ),
        "image_uri": "docker://tensorflow/tensorflow:2.8.2-gpu",
        "installed_chrombpnet": False,
        "runtime_id": "chrombpnet_native_tf2_8_2",
        "tensorflow_version": probe["tensorflow_version"],
    },
)
PY

mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet native runtime acquisition: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
