#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-overlay
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${RUNTIME}/tensorflow-tensorflow-2.8.2-gpu.sif"
RUNTIME_ARTIFACTS_SHA256=b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229
IMAGE_SHA256=96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa
TARGET="${ROOT}/executions/environments/chrombpnet-native-core-overlay-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite runtime overlay: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/site-packages" \
  "${STAGE}/source-distributions" \
  "${STAGE}/wheels" \
  "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${IMAGE_SHA256}" "${IMAGE}" | sha256sum --check --strict

cat > "${STAGE}/requirements.lock" <<'EOF'
absl-py==1.0.0
cloudpickle==2.2.1
decorator==5.1.1
dm-tree==0.1.8
gast==0.5.3
joblib==1.2.0
numpy==1.23.4
pandas==1.3.5
protobuf==3.20.0
pyBigWig==0.3.22
python-dateutil==2.8.2
pytz==2022.1
scikit-learn==1.1.3
scipy==1.8.1
six==1.14.0
tensorflow-probability==0.15.0
threadpoolctl==3.1.0
tqdm==4.48.2
EOF

module load singularity/4.2.2
module load python/3.11.5-GCCcore-13.2.0
singularity exec --cleanenv --bind "${STAGE}:${STAGE}:rw" "${IMAGE}" \
  python -m pip download \
    --disable-pip-version-check \
    --only-binary=:all: \
    --dest "${STAGE}/wheels" \
    --requirement "${STAGE}/requirements.lock"
singularity exec --cleanenv --bind "${STAGE}:${STAGE}:rw" "${IMAGE}" \
  python -m pip download \
    --disable-pip-version-check \
    --no-binary=pyfaidx \
    --no-deps \
    --dest "${STAGE}/source-distributions" \
    pyfaidx==0.6.1
singularity exec --cleanenv --bind "${STAGE}:${STAGE}:rw" "${IMAGE}" \
  python -m pip wheel \
    --disable-pip-version-check \
    --no-deps \
    --wheel-dir "${STAGE}/wheels" \
    "${STAGE}"/source-distributions/pyfaidx-0.6.1.tar.gz
singularity exec --cleanenv --bind "${STAGE}:${STAGE}:rw" "${IMAGE}" \
  python -m pip install \
    --disable-pip-version-check \
    --no-index \
    --no-deps \
    --target "${STAGE}/site-packages" \
    "${STAGE}"/wheels/*.whl

sha256sum "${STAGE}"/wheels/*.whl > "${STAGE}/validation/wheels.sha256"
sha256sum "${STAGE}"/source-distributions/* \
  > "${STAGE}/validation/source-distributions.sha256"
singularity exec --cleanenv \
  --bind "${STAGE}:${STAGE}:ro" \
  --env "PYTHONPATH=${STAGE}/site-packages" \
  "${IMAGE}" python - <<'PY' > "${STAGE}/validation/import_probe.json"
import json
from importlib.metadata import version
import sys

import cloudpickle
import google.protobuf
import joblib
import numpy
import pandas
import pyBigWig
import pyfaidx
import scipy
import sklearn
import tensorflow as tf
import tensorflow_probability as tfp

expected = {
    "absl_py": "1.0.0",
    "cloudpickle": "2.2.1",
    "gast": "0.5.3",
    "joblib": "1.2.0",
    "numpy": "1.23.4",
    "pandas": "1.3.5",
    "protobuf": "3.20.0",
    "pyBigWig": "0.3.22",
    "pyfaidx": "0.6.1",
    "scipy": "1.8.1",
    "scikit_learn": "1.1.3",
    "six": "1.14.0",
    "tensorflow": "2.8.2",
    "tensorflow_probability": "0.15.0",
}
observed = {
    "absl_py": version("absl-py"),
    "cloudpickle": cloudpickle.__version__,
    "gast": version("gast"),
    "joblib": joblib.__version__,
    "numpy": numpy.__version__,
    "pandas": pandas.__version__,
    "protobuf": google.protobuf.__version__,
    "pyBigWig": pyBigWig.__version__,
    "pyfaidx": pyfaidx.__version__,
    "scipy": scipy.__version__,
    "scikit_learn": sklearn.__version__,
    "six": version("six"),
    "tensorflow": tf.__version__,
    "tensorflow_probability": tfp.__version__,
}
if observed != expected:
    raise SystemExit(f"runtime package versions differ: {observed!r}")
print(
    json.dumps(
        {
            "python": sys.version,
            "versions": observed,
            "gpu_visible_during_cpu_build": bool(
                tf.config.list_physical_devices("GPU")
            ),
        },
        sort_keys=True,
        separators=(",", ":"),
    )
)
PY
singularity exec --cleanenv \
  --bind "${STAGE}:${STAGE}:ro" \
  --env "PYTHONPATH=${STAGE}/site-packages" \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
singularity --version > "${STAGE}/validation/singularity.version.txt"

cat > "${STAGE}/SOURCE.tsv" <<EOF
field\tvalue
runtime_id\tchrombpnet_native_core_tf2_8_2
base_runtime_artifacts_sha256\t${RUNTIME_ARTIFACTS_SHA256}
base_image_sha256\t${IMAGE_SHA256}
python\t3.8.10
tensorflow\t2.8.2
tensorflow_probability\t0.15.0
dependency_policy\tpinned_core_training_dependencies_only
EOF

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "chrombpnet_native_core_overlay",
        "base_image_sha256": "96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa",
        "base_runtime_artifacts_sha256": "b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229",
        "requirements_lock_sha256": sha256_file(root / "requirements.lock"),
        "runtime_id": "chrombpnet_native_core_tf2_8_2",
    },
)
PY

mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet native core overlay: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
