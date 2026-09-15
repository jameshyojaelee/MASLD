#!/usr/bin/env bash
#SBATCH --job-name=model-runtime-142
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/environments/scvi-tools-1.4.2-inductive-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
EXPECTED_SCVI_WHEEL_SHA256=c453b75b0fa0d222bf0eb3a531607d17f5ea0dbe5b848ab12fd8cc24599884e5

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite scvi-tools 1.4.2 runtime: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/downloads" \
  "${STAGE}/validation" \
  "${STAGE}/validation/xdg-cache"

module purge
module load python/3.12.3-GCCcore-13.3.0
python -m venv "${STAGE}/venv"
PYTHON=${STAGE}/venv/bin/python
export PIP_NO_CACHE_DIR=1
export XDG_CACHE_HOME=${STAGE}/validation/xdg-cache
"${PYTHON}" -m pip install \
  --disable-pip-version-check \
  --upgrade \
  pip==25.2 \
  setuptools==80.9.0 \
  wheel==0.45.1
"${PYTHON}" -m pip download \
  --disable-pip-version-check \
  --no-deps \
  --only-binary=:all: \
  --dest "${STAGE}/downloads" \
  scvi-tools==1.4.2
SCVI_WHEELS=("${STAGE}"/downloads/scvi_tools-1.4.2-*.whl)
if [[ "${#SCVI_WHEELS[@]}" -ne 1 || ! -f "${SCVI_WHEELS[0]}" ]]; then
  printf 'Expected exactly one scvi-tools 1.4.2 wheel.\n' >&2
  exit 1
fi
printf '%s  %s\n' "${EXPECTED_SCVI_WHEEL_SHA256}" "${SCVI_WHEELS[0]}" \
  | sha256sum --check --strict
"${PYTHON}" -m pip install \
  --disable-pip-version-check \
  --report "${STAGE}/validation/pip-install-report.json" \
  "${SCVI_WHEELS[0]}"
"${PYTHON}" -m pip check > "${STAGE}/validation/pip-check.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" -m pip list --format=json > "${STAGE}/validation/pip-list.json"

cd "${ROOT}"
export PYTHONPATH=${ROOT}/src
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${PYTHON}" -m unittest tests.unit.test_scvi_scanvi_inductive \
  > "${STAGE}/validation/unit-tests.txt" \
  2>&1
"${PYTHON}" scripts/probe_scvi_scanvi_1_4_2_inductive.py \
  --output "${STAGE}/validation/inductive-probe" \
  > "${STAGE}/validation/inductive-probe.stdout.json" \
  2> "${STAGE}/validation/inductive-probe.stderr.txt"
sha256sum \
  slurm/acquire_scvi_tools_1_4_2_inductive_runtime.sbatch \
  scripts/probe_scvi_scanvi_1_4_2_inductive.py \
  src/masld_bench/adapters/scvi_scanvi_inductive.py \
  tests/unit/test_scvi_scanvi_inductive.py \
  "${SCVI_WHEELS[0]}" \
  > "${STAGE}/source.sha256"

PROBE=${STAGE}/validation/inductive-probe/probe.json
"${PYTHON}" - "${PROBE}" "${EXPECTED_SCVI_WHEEL_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

probe = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
if (
    probe.get("status") != "pass"
    or probe.get("scvi_tools") != "1.4.2"
    or probe.get("query_training_executed") is not False
    or probe.get("query_adaptation_executed") is not False
    or probe.get("direct_scvi_inference") is not True
    or probe.get("direct_scanvi_inference") is not True
    or probe.get("scvi_state_unchanged") is not True
    or probe.get("scanvi_state_unchanged") is not True
    or probe.get("project_data_read") is not False
    or probe.get("sealed_outcomes_read") is not False
):
    raise SystemExit("scvi-tools 1.4.2 inductive admission probe differs")
print(json.dumps({"probe": probe, "wheel_sha256": sys.argv[2]}, sort_keys=True))
PY

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${STAGE}" "${EXPECTED_SCVI_WHEEL_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads(
    (root / "validation/inductive-probe/probe.json").read_text(encoding="utf-8")
)
freeze_tree(
    root / "validation/inductive-probe",
    {
        "artifact_class": "scvi_scanvi_inductive_query_probe",
        "scvi_tools_version": "1.4.2",
        "query_training_executed": False,
        "query_adaptation_executed": False,
        "learned_state_unchanged": True,
        "project_data_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
freeze_tree(
    root,
    {
        "artifact_class": "scvi_tools_1_4_2_inductive_runtime",
        "model_roster": ["scvi_baseline", "scanvi_baseline"],
        "scvi_tools_version": "1.4.2",
        "scvi_tools_wheel_sha256": sys.argv[2],
        "code_commit": "dc95497",
        "code_license": "BSD-3-Clause",
        "python_version": probe["python"],
        "direct_inductive_query_probe": "passed",
        "query_training_allowed": False,
        "query_adaptation_allowed": False,
        "project_data_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'scvi-tools 1.4.2 inductive runtime: %s\n' "${TARGET}"
sha256sum \
  "${TARGET}/ARTIFACTS.json" \
  "${TARGET}/validation/inductive-probe/ARTIFACTS.json"
