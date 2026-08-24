#!/usr/bin/env bash
#SBATCH --job-name=masld-midas030-runtime
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
PACKAGE=scmidas
VERSION=0.3.0
WHEEL_SHA256=24df81372decc6d3edf6c2481d7b0c58f4685386d2ac938564946916f1bd9522
TARGET="${ROOT}/executions/environments/midas-v0.3.0-portable-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed MIDAS runtime: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite MIDAS runtime: %s\n' "${TARGET}" >&2
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
WHEEL=$(find "${STAGE}/source" -maxdepth 1 -type f -name 'scmidas-0.3.0-*.whl' -print)
if [[ -z "${WHEEL}" || $(printf '%s\n' "${WHEEL}" | wc -l) -ne 1 ]]; then
  printf 'MIDAS wheel census differs\n' >&2
  exit 1
fi
printf '%s  %s\n' "${WHEEL_SHA256}" "${WHEEL}" | sha256sum --check --strict
PYTHONPATH="${TORCH}/site-packages:${BASE}/site-packages" python -m pip install \
  --disable-pip-version-check --no-deps --target "${STAGE}/site-packages" "${WHEEL}" \
  > "${STAGE}/validation/install.stdout.txt" \
  2> "${STAGE}/validation/install.stderr.txt"
PYTHONPATH="${TORCH}/site-packages:${BASE}/site-packages" python -m pip download \
  --disable-pip-version-check --no-deps --only-binary=:all: \
  --dest "${STAGE}/source" "toml==0.10.2" \
  > "${STAGE}/validation/dependency-download.stdout.txt" \
  2> "${STAGE}/validation/dependency-download.stderr.txt"
TOML_WHEEL=$(find "${STAGE}/source" -maxdepth 1 -type f -name 'toml-0.10.2-*.whl' -print)
if [[ -z "${TOML_WHEEL}" || $(printf '%s\n' "${TOML_WHEEL}" | wc -l) -ne 1 ]]; then
  printf 'MIDAS toml wheel census differs\n' >&2
  exit 1
fi
PYTHONPATH="${TORCH}/site-packages:${BASE}/site-packages" python -m pip install \
  --disable-pip-version-check --no-deps --target "${STAGE}/site-packages" \
  "${TOML_WHEEL}" \
  > "${STAGE}/validation/dependency-install.stdout.txt" \
  2> "${STAGE}/validation/dependency-install.stderr.txt"
PYTHONPATH="${STAGE}/site-packages:${TORCH}/site-packages:${BASE}/site-packages" \
  python - "${STAGE}" <<'PY' > "${STAGE}/validation/import.json"
from importlib.metadata import version
import json
from pathlib import Path
import sys

root = Path(sys.argv[1]).resolve(strict=True)
import scmidas
from scmidas.model import MIDAS

receipt = {
    "schema_version": "masld-bench-midas-runtime-v1",
    "status": "pass",
    "distribution": "scmidas",
    "version": version("scmidas"),
    "midas_class_imported": MIDAS.__name__ == "MIDAS",
    "project_data_read": False,
    "outcomes_read": False,
    "model_fit_executed": False,
    "unsafe_checkpoint_load_executed": False,
}
if receipt["version"] != "0.3.0" or not receipt["midas_class_imported"]:
    raise SystemExit("MIDAS runtime contract differs")
print(json.dumps(receipt, sort_keys=True))
(root / "validation/receipt.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
PY
PYTHONPATH="${STAGE}/site-packages:${TORCH}/site-packages:${BASE}/site-packages" \
  python -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum "${ROOT}/slurm/midas_prepare_v0_3_0_runtime.sbatch" > "${STAGE}/source.sha256"
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
    or receipt["version"] != "0.3.0"
    or not receipt["midas_class_imported"]
    or receipt["project_data_read"]
    or receipt["outcomes_read"]
    or receipt["model_fit_executed"]
    or receipt["unsafe_checkpoint_load_executed"]
):
    raise SystemExit("MIDAS runtime receipt differs")
freeze_tree(root, {
    "artifact_class": "midas_v0_3_0_portable_runtime",
    "distribution": "scmidas",
    "version": "0.3.0",
    "wheel_sha256": "24df81372decc6d3edf6c2481d7b0c58f4685386d2ac938564946916f1bd9522",
    "project_data_read": False,
    "outcomes_read": False,
    "model_fit_executed": False,
    "unsafe_checkpoint_load_executed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'MIDAS v0.3.0 portable runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
