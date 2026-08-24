#!/usr/bin/env bash
#SBATCH --job-name=masld-scbutterfly009-runtime
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
BASE_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_ARTIFACTS_SHA256=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
PACKAGE=scButterfly
VERSION=0.0.9
WHEEL_SHA256=98aab0c1d8445f69a08126821d6a3d250d17c3d2065a35b339e711754f1293e7
TARGET="${ROOT}/executions/environments/scbutterfly-v0.0.9-portable-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scButterfly runtime: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scButterfly runtime: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/site-packages" "${STAGE}/source" "${STAGE}/validation"
for record in \
  "${BASE_ARTIFACTS_SHA256} ${BASE}/ARTIFACTS.json" \
  "${TORCH_ARTIFACTS_SHA256} ${TORCH}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BASE}" "${TORCH}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${TORCH}/site-packages:${BASE}/site-packages" python -m pip download \
  --disable-pip-version-check --no-deps --only-binary=:all: \
  --dest "${STAGE}/source" "${PACKAGE}==${VERSION}" \
  > "${STAGE}/validation/download.stdout.txt" \
  2> "${STAGE}/validation/download.stderr.txt"
WHEEL=$(find "${STAGE}/source" -maxdepth 1 -type f -iname 'scbutterfly-0.0.9-*.whl' -print)
if [[ -z "${WHEEL}" || $(printf '%s\n' "${WHEEL}" | wc -l) -ne 1 ]]; then
  printf 'scButterfly wheel census differs\n' >&2
  exit 1
fi
printf '%s  %s\n' "${WHEEL_SHA256}" "${WHEEL}" | sha256sum --check --strict
PYTHONPATH="${TORCH}/site-packages:${BASE}/site-packages" python -m pip install \
  --disable-pip-version-check --no-deps --target "${STAGE}/site-packages" "${WHEEL}" \
  > "${STAGE}/validation/install.stdout.txt" \
  2> "${STAGE}/validation/install.stderr.txt"
PYTHONPATH="${STAGE}/site-packages:${TORCH}/site-packages:${BASE}/site-packages" \
  python - "${STAGE}" <<'PY' > "${STAGE}/validation/import.json"
from importlib.metadata import version
import importlib.util
import json
from pathlib import Path
import sys

root = Path(sys.argv[1]).resolve(strict=True)
component_path = root / "site-packages/scButterfly/model_component.py"
spec = importlib.util.spec_from_file_location("scbutterfly_model_component", component_path)
if spec is None or spec.loader is None:
    raise SystemExit("scButterfly model component spec differs")
component = importlib.util.module_from_spec(spec)
spec.loader.exec_module(component)

receipt = {
    "schema_version": "masld-bench-scbutterfly-runtime-v1",
    "status": "pass",
    "distribution": "scButterfly",
    "version": version("scButterfly"),
    "translator_class_imported": component.Translator.__name__ == "Translator",
    "high_level_leaky_api_imported": False,
    "project_data_read": False,
    "outcomes_read": False,
    "model_fit_executed": False,
}
if receipt["version"] != "0.0.9" or not receipt["translator_class_imported"]:
    raise SystemExit("scButterfly runtime contract differs")
print(json.dumps(receipt, sort_keys=True))
(root / "validation/receipt.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
PY
PYTHONPATH="${STAGE}/site-packages:${TORCH}/site-packages:${BASE}/site-packages" \
  python -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum "${ROOT}/slurm/scbutterfly_prepare_v0_0_9_runtime.sbatch" > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["version"] != "0.0.9"
    or not receipt["translator_class_imported"]
    or receipt["high_level_leaky_api_imported"]
    or receipt["project_data_read"]
    or receipt["outcomes_read"]
    or receipt["model_fit_executed"]
):
    raise SystemExit("scButterfly runtime receipt differs")
freeze_tree(root, {
    "artifact_class": "scbutterfly_v0_0_9_portable_runtime",
    "distribution": "scButterfly",
    "version": "0.0.9",
    "wheel_sha256": "98aab0c1d8445f69a08126821d6a3d250d17c3d2065a35b339e711754f1293e7",
    "project_data_read": False,
    "outcomes_read": False,
    "model_fit_executed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scButterfly v0.0.9 portable runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
