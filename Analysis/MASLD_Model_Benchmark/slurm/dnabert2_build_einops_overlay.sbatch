#!/usr/bin/env bash
#SBATCH --job-name=masld-dnabert2-einops
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
TARGET="${ROOT}/executions/environments/dnabert2-einops-0.6.1-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite DNABERT-2 overlay: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/site-packages" "${STAGE}/wheels" "${STAGE}/validation"
"${BASE}" -m pip download \
  --disable-pip-version-check \
  --only-binary=:all: \
  --no-deps \
  --dest "${STAGE}/wheels" \
  einops==0.6.1
"${BASE}" -m pip install \
  --disable-pip-version-check \
  --no-index \
  --no-deps \
  --target "${STAGE}/site-packages" \
  "${STAGE}"/wheels/*.whl
sha256sum "${STAGE}"/wheels/*.whl > "${STAGE}/validation/wheels.sha256"
PYTHONPATH="${STAGE}/site-packages" "${BASE}" - <<'PY' \
  > "${STAGE}/validation/import_probe.json"
import json
from importlib.metadata import distribution

package = distribution("einops")
if package.version != "0.6.1":
    raise SystemExit("einops version differs")
print(
    json.dumps(
        {
            "name": "einops",
            "version": package.version,
            "license": package.metadata.get("License"),
            "custom_model_code_imported": False,
            "checkpoint_loaded": False,
        },
        sort_keys=True,
    )
)
PY
PYTHONPATH="${STAGE}/site-packages" "${BASE}" -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum "${ROOT}/slurm/dnabert2_build_einops_overlay.sbatch" \
  > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation" / "import_probe.json").read_text())
if probe["version"] != "0.6.1" or probe["checkpoint_loaded"]:
    raise SystemExit("DNABERT-2 overlay receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "dnabert2_einops_overlay",
        "einops": "0.6.1",
        "wheel_sha256": sha256_file(next((root / "wheels").glob("*.whl"))),
        "checkpoint_loaded": False,
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'DNABERT-2 einops overlay: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
