#!/usr/bin/env bash
#SBATCH --job-name=model-check-328
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/model-check-328-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
RUNTIME="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738"
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON="${RUNTIME}/env/bin/python"

preserve_failure() {
  status=$?
  if [[ -d "${STAGE}" ]]; then
    if [[ "${status}" -ne 0 && ! -e "${TARGET}.failed" ]]; then
      mv "${STAGE}" "${TARGET}.failed"
    else
      rm -rf -- "${STAGE}"
    fi
  fi
  return "${status}"
}
trap preserve_failure EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite current Corgi+ readiness audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/source-cache/corgi-current/corgi" \
  "${STAGE}/source-cache/reproduction-current/code" \
  "${STAGE}/validation" "${STAGE}/runtime"

printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
curl_flags=(--fail --location --silent --show-error --retry 6 --retry-all-errors --connect-timeout 30)

curl "${curl_flags[@]}" --output "${STAGE}/source-cache/corgi.repo.json" \
  https://api.github.com/repos/ekinda/corgi
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/reproduction.repo.json" \
  https://api.github.com/repos/ekinda/corgi-reproduction
read -r corgi_branch reproduction_branch < <(
  "${PYTHON}" - "${STAGE}/source-cache/corgi.repo.json" \
    "${STAGE}/source-cache/reproduction.repo.json" <<'PY'
from pathlib import Path
import json
import sys

values = [json.loads(Path(path).read_text())["default_branch"] for path in sys.argv[1:]]
print(*values)
PY
)
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/corgi.commit.json" \
  "https://api.github.com/repos/ekinda/corgi/commits/${corgi_branch}"
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/reproduction.commit.json" \
  "https://api.github.com/repos/ekinda/corgi-reproduction/commits/${reproduction_branch}"
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/corgi.tags.json" \
  https://api.github.com/repos/ekinda/corgi/tags?per_page=100
read -r corgi_revision reproduction_revision < <(
  "${PYTHON}" - "${STAGE}/source-cache/corgi.commit.json" \
    "${STAGE}/source-cache/reproduction.commit.json" <<'PY'
from pathlib import Path
import json
import sys

values = [json.loads(Path(path).read_text())["sha"] for path in sys.argv[1:]]
print(*values)
PY
)

for relative in LICENSE pyproject.toml corgi/predict.py corgi/config_corgiplus.py corgi/trainer_corgiplus.py; do
  curl "${curl_flags[@]}" \
    --output "${STAGE}/source-cache/corgi-current/${relative}" \
    "https://raw.githubusercontent.com/ekinda/corgi/${corgi_revision}/${relative}"
done
curl "${curl_flags[@]}" \
  --output "${STAGE}/source-cache/reproduction-current/code/corgiplus_delta_bench_tta.py" \
  "https://raw.githubusercontent.com/ekinda/corgi-reproduction/${reproduction_revision}/code/corgiplus_delta_bench_tta.py"
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/zenodo-18630048.json" \
  https://zenodo.org/api/records/18630048
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/zenodo-20686797.json" \
  https://zenodo.org/api/records/20686797

export PYTHONPATH="${ROOT}/src"
cd "${ROOT}"
"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_corgi_plus_current_terminal_readiness \
  > "${STAGE}/validation/unit-test.stdout.txt" 2>&1
"${PYTHON}" "${ROOT}/scripts/audit_corgi_plus_current_terminal_readiness.py" \
  --project-root "${ROOT}" \
  --source-cache "${STAGE}/source-cache" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json"

"${PYTHON}" - "${STAGE}/audit/audit.json" <<'PY'
from pathlib import Path
import json
import sys

audit = json.loads(Path(sys.argv[1]).read_text())
assert audit["audit_status"] == "passed_terminal_blocked"
assert audit["checkpoint"]["sha256"] == "01855614bfeffa72d5bcaaa28729c022ccfbc8f5fdcc0add8e04ba4533e94ccd"
assert audit["checkpoint"]["auxiliary_width"] == 26
assert audit["licenses"]["code_conflict"] == "UNRESOLVED"
assert not audit["licenses"]["open_champion_eligible"]
assert not audit["native_input_contract"]["native_input_bundle_released"]
assert not audit["native_input_contract"]["gene_count_matrix_is_valid_local_RNA"]
assert audit["reference_contract"]["biological_execution_gate"] == "blocked"
assert not any(row["corgi_plus_family_native_input_complete"] for row in audit["project_fixture_fit"].values())
assert audit["nearest_executable_lane"]["model_id"] == "corgi_regular"
assert not audit["nearest_executable_lane"]["relabel_as_Corgi_plus"]
assert not audit["context_borzoi_boundary"]["nearest_executable_now"]
assert not audit["terminal_disposition"]["production_or_GPU_item_authorized"]
assert not audit["checkpoint_deserialized_by_this_audit"]
assert not audit["model_forward_executed_by_this_audit"]
assert not audit["biological_data_read"]
assert not audit["development_outcomes_read"]
assert not audit["sealed_features_read"]
assert not audit["sealed_outcomes_read"]
assert not audit["gpu_job_submitted"]
PY

"${PYTHON}" --version > "${STAGE}/runtime/python-version.txt" 2>&1
"${PYTHON}" -m pip freeze --all > "${STAGE}/runtime/pip-freeze.txt"
uname -a > "${STAGE}/runtime/uname.txt"
env | LC_ALL=C sort | grep '^SLURM_' > "${STAGE}/runtime/slurm-environment.txt"
sha256sum \
  "${ROOT}/scripts/audit_corgi_plus_current_terminal_readiness.py" \
  "${ROOT}/tests/unit/test_audit_corgi_plus_current_terminal_readiness.py" \
  "${ROOT}/slurm/audit_corgi_plus_current_terminal_readiness_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "audit/audit.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "corgi_plus_current_terminal_readiness",
        "model": "corgi_plus",
        "audit_status": audit["audit_status"],
        "released_biological_execution_authorized": False,
        "nearest_executable_lane": "corgi_regular",
        "context_borzoi_build_authorized": False,
        "checkpoint_deserialized": False,
        "model_forward_executed": False,
        "biological_data_read": False,
        "development_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_outcomes_read": False,
        "gpu_job_submitted": False,
        "status": "passed_terminal_blocked",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Current Corgi+ terminal readiness: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
