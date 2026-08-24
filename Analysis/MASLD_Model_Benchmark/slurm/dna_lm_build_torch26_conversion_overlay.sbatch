#!/usr/bin/env bash
#SBATCH --job-name=masld-evo-convert-runtime
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
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
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
BASE_SHA=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_SHA=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
TARGET="${ROOT}/executions/environments/dna-lm-torch2.6-conversion-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite torch2.6 conversion overlay: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/site-packages" "${STAGE}/wheels" "${STAGE}/validation"
printf '%s  %s\n' "${BASE_SHA}" "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${TORCH_SHA}" "${TORCH}/ARTIFACTS.json" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BASE}" "${TORCH}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
RUNTIME_PATH="${STAGE}/site-packages:${TORCH}/site-packages:${BASE}/site-packages"
PYTHONPATH="${TORCH}/site-packages:${BASE}/site-packages" python -m pip download \
  --disable-pip-version-check \
  --only-binary=:all: \
  --no-deps \
  --dest "${STAGE}/wheels" \
  safetensors==0.5.2
PYTHONPATH="${TORCH}/site-packages:${BASE}/site-packages" python -m pip install \
  --disable-pip-version-check \
  --no-index \
  --no-deps \
  --target "${STAGE}/site-packages" \
  "${STAGE}"/wheels/*.whl
sha256sum "${STAGE}"/wheels/*.whl > "${STAGE}/validation/wheels.sha256"
PYTHONPATH="${RUNTIME_PATH}" python - <<'PY' > "${STAGE}/validation/import_probe.json"
import json
from importlib.metadata import distribution
import sys
import torch
import safetensors

package = distribution("safetensors")
if (
    sys.version_info[:2] != (3, 12)
    or torch.__version__ != "2.6.0+cu124"
    or package.version != "0.5.2"
    or not hasattr(torch.serialization, "add_safe_globals")
):
    raise SystemExit("torch2.6 conversion runtime differs")
print(
    json.dumps(
        {
            "python": sys.version,
            "torch": torch.__version__,
            "torch_cuda": torch.version.cuda,
            "safetensors": package.version,
            "weights_only_safe_globals_api": True,
            "checkpoint_loaded": False,
            "model_forward_executed": False,
        },
        sort_keys=True,
    )
)
PY
PYTHONPATH="${RUNTIME_PATH}" python -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/slurm/dna_lm_build_torch26_conversion_overlay.sbatch" \
  "${BASE}/ARTIFACTS.json" \
  "${TORCH}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation" / "import_probe.json").read_text())
if probe["checkpoint_loaded"] or not probe["weights_only_safe_globals_api"]:
    raise SystemExit("torch2.6 conversion runtime receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "dna_lm_torch2_6_conversion_overlay",
        "base_artifacts_sha256": "30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af",
        "torch_artifacts_sha256": "d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3",
        "torch": "2.6.0+cu124",
        "safetensors": "0.5.2",
        "wheel_sha256": sha256_file(next((root / "wheels").glob("*.whl"))),
        "checkpoint_loaded": False,
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'DNA-LM torch2.6 conversion overlay: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
