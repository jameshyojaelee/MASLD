#!/usr/bin/env bash
#SBATCH --job-name=model-data-ep901
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python
TARGET="${ROOT}/executions/model-data-ep901-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
export PYTHONPATH="${ROOT}:${ROOT}/src"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite EpiBERT primary evidence: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

"${PYTHON}" -m unittest -v \
  tests.unit.test_acquire_epibert_primary_evidence \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/acquire_epibert_primary_evidence.py" \
  --repository-url https://github.com/naumanjaved/EpiBERT.git \
  --pinned-revision dae61b434f885991c29e2fde3fbf7f70d0ff9b80 \
  --tag v1.0.0 \
  --pmc-id PMC11872434 \
  --pdftotext /gpfs/commons/home/jameslee/.local/bin/pdftotext \
  --output "${STAGE}/sources" \
  > "${STAGE}/validation/acquisition.stdout.json" \
  2> "${STAGE}/validation/acquisition.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
/usr/bin/git --version > "${STAGE}/validation/git.version.txt"
/usr/bin/curl --version > "${STAGE}/validation/curl.version.txt"
/gpfs/commons/home/jameslee/.local/bin/pdftotext -v \
  > "${STAGE}/validation/pdftotext.version.txt" 2>&1
sha256sum \
  "${ROOT}/scripts/acquire_epibert_primary_evidence.py" \
  "${ROOT}/tests/unit/test_acquire_epibert_primary_evidence.py" \
  "${ROOT}/slurm/acquire_epibert_primary_evidence_io.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${TARGET}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace


stage = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2])
receipt = json.loads((stage / "sources/acquisition_receipt.json").read_text(encoding="utf-8"))
tests = (stage / "validation/tests.stderr.txt").read_text(encoding="utf-8")
if (
    receipt["status"] != "pass"
    or receipt["repository"]["pinned_revision"] != "dae61b434f885991c29e2fde3fbf7f70d0ff9b80"
    or receipt["repository"]["tag"] != "v1.0.0"
    or receipt["paper"]["pmc_id"] != "PMC11872434"
    or receipt["project_data_read"]
    or receipt["outcomes_read"]
    or receipt["sealed_assets_read"]
    or "Ran 3 tests" not in tests
    or "OK" not in tests
):
    raise SystemExit("EpiBERT primary-evidence acquisition differs")
freeze_tree(
    stage,
    {
        "artifact_class": "epibert_primary_source_history_and_paper_evidence",
        "model_id": "epibert",
        "pinned_revision": "dae61b434f885991c29e2fde3fbf7f70d0ff9b80",
        "tag": "v1.0.0",
        "pmc_id": "PMC11872434",
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "status": "pass",
    },
)
publish_directory_noreplace(stage, target)
PY

trap - EXIT
printf 'EpiBERT primary evidence: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
