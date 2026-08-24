#!/usr/bin/env bash
#SBATCH --job-name=masld-epibert-tf213-runtime
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/environments/epibert-tensorflow-2.13.0-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
IMAGE="${STAGE}/tensorflow-tensorflow-2.13.0-gpu.sif"
OCI_DIGEST=sha256:b4676741c491bff3d0f29c38c369281792c7d5c5bfa2b1aa93e5231a8d236323

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite EpiBERT runtime artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

module load singularity/4.2.2
singularity pull --disable-cache "${IMAGE}" \
  "docker://tensorflow/tensorflow@${OCI_DIGEST}"
singularity inspect --json "${IMAGE}" > "${STAGE}/validation/singularity-inspect.json"
singularity exec --cleanenv "${IMAGE}" python - <<'PY' \
  > "${STAGE}/validation/runtime-probe.json"
import json
import platform
import tensorflow as tf

if tf.__version__ != "2.13.0":
    raise SystemExit(f"TensorFlow version differs: {tf.__version__}")
print(json.dumps({
    "python": platform.python_version(),
    "tensorflow": tf.__version__,
    "gpu_visible_during_cpu_acquisition": bool(tf.config.list_physical_devices("GPU")),
    "status": "pass",
}, sort_keys=True))
PY
sha256sum "${IMAGE}" > "${STAGE}/validation/image.sha256"
sha256sum "${ROOT}/slurm/epibert_acquire_tf213_runtime.sbatch" > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${OCI_DIGEST}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
image = root / "tensorflow-tensorflow-2.13.0-gpu.sif"
receipt = json.loads((root / "validation/runtime-probe.json").read_text())
if receipt.get("status") != "pass" or receipt.get("tensorflow") != "2.13.0":
    raise SystemExit("EpiBERT TensorFlow runtime probe differs")
freeze_tree(
    root,
    {
        "artifact_class": "epibert_tensorflow_runtime",
        "container_source": f"docker://tensorflow/tensorflow@{sys.argv[2]}",
        "oci_manifest_digest": sys.argv[2],
        "image_sha256": sha256_file(image),
        "tensorflow": "2.13.0",
        "tensorflow_addons_installed": False,
        "full_epibert_architecture_import_ready": False,
        "no_outcome_tensor_fixture_ready": True,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'EpiBERT TensorFlow 2.13 runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
