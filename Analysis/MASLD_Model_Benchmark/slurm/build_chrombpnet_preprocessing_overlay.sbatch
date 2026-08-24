#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-prep-overlay
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
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
LOCK="${ROOT}/config/runtimes/chrombpnet_preprocessing_requirements.lock"
TARGET="${ROOT}/executions/environments/chrombpnet-native-preprocessing-overlay-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite preprocessing overlay: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/site-packages" "${STAGE}/wheels" \
  "${STAGE}/validation" "${STAGE}/mplconfig"
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
cp "${LOCK}" "${STAGE}/requirements.lock"

module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${STAGE}:${STAGE}:rw" \
  "${IMAGE}" python -m pip download \
    --disable-pip-version-check \
    --only-binary=:all: \
    --no-deps \
    --dest "${STAGE}/wheels" \
    --requirement "${STAGE}/requirements.lock"
singularity exec --cleanenv \
  --bind "${STAGE}:${STAGE}:rw" \
  "${IMAGE}" python -m pip install \
    --disable-pip-version-check \
    --no-index \
    --no-deps \
    --target "${STAGE}/site-packages" \
    "${STAGE}"/wheels/*.whl
sha256sum "${STAGE}"/wheels/*.whl > "${STAGE}/validation/wheels.sha256"

singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${STAGE}/site-packages:${CORE}/site-packages" \
  --env "MPLCONFIGDIR=${STAGE}/mplconfig" \
  --env MPLBACKEND=Agg \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from importlib.metadata import metadata, version
from pathlib import Path
import json
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy
import pandas

root = Path(sys.argv[1]).resolve(strict=True)
expected = {
    "cycler": "0.11.0",
    "fonttools": "4.38.0",
    "kiwisolver": "1.4.4",
    "matplotlib": "3.5.3",
    "packaging": "21.3",
    "Pillow": "9.3.0",
    "pyparsing": "3.0.9",
}
observed = {name: version(name) for name in expected}
if observed != expected:
    raise SystemExit(f"preprocessing package versions differ: {observed!r}")
if numpy.__version__ != "1.23.4" or pandas.__version__ != "1.3.5":
    raise SystemExit("core numerical runtime changed under preprocessing overlay")
figure = root / "validation" / "matplotlib-probe.png"
plt.hist([[0.4, 0.5], [0.5, 0.6]], bins=5, density=True)
plt.savefig(figure)
plt.close()
if not figure.is_file() or figure.stat().st_size == 0:
    raise SystemExit("matplotlib Agg probe did not emit a figure")
licenses = {name: metadata(name).get("License") for name in expected}
result = {
    "schema_version": "masld-bench-chrombpnet-preprocessing-overlay-v1",
    "status": "pass",
    "versions": observed,
    "licenses": licenses,
    "numpy": numpy.__version__,
    "pandas": pandas.__version__,
    "matplotlib_backend": matplotlib.get_backend(),
    "figure_size_bytes": figure.stat().st_size,
}
(root / "validation" / "import-probe.json").write_text(
    json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(result, sort_keys=True))
PY
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:ro" \
  --env "PYTHONPATH=${STAGE}/site-packages:${CORE}/site-packages" \
  --env "MPLCONFIGDIR=${STAGE}/mplconfig" \
  --env MPLBACKEND=Agg \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
singularity --version > "${STAGE}/validation/singularity.version.txt"
sha256sum \
  "${ROOT}/slurm/build_chrombpnet_preprocessing_overlay.sbatch" \
  "${LOCK}" > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "chrombpnet_native_preprocessing_overlay",
            "runtime_id": "chrombpnet_native_preprocessing_tf2_8_2",
            "base_image_sha256": "96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa",
            "core_overlay_artifacts_sha256": "61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1",
            "requirements_lock_sha256": sha256_file(root / "requirements.lock"),
            "matplotlib": "3.5.3",
            "purpose": "chrombpnet_native_preprocessing_only",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet preprocessing overlay: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
