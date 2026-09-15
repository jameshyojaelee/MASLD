#!/usr/bin/env bash
#SBATCH --job-name=model-check-620
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SCRIPT="${ROOT}/scripts/aggregate_lsgkm_gse281364_dinucleotide_scores.py"
TEST="${ROOT}/tests/unit/test_aggregate_lsgkm_gse281364_dinucleotide_scores.py"
TARGET="${ROOT}/executions/model-check-620-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite aggregation tests: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' 3a0367a350893989228536478acd83dcd052063342b5f9fc5fb67195452efc89 "${SCRIPT}" | sha256sum --check --strict
printf '%s  %s\n' 7f9e0b9365409f16202b7de27ab892cd9bc0bedbd4f2f0183992b2543b047337 "${TEST}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}" python -m unittest -v \
  tests.unit.test_aggregate_lsgkm_gse281364_dinucleotide_scores \
  > "${STAGE}/tests.stdout.txt" 2> "${STAGE}/tests.stderr.txt"
python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
sha256sum "${SCRIPT}" "${TEST}" \
  "${ROOT}/slurm/test_lsgkm_gse281364_aggregation_cpu.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = {
    "schema_version": "masld-bench-lsgkm-aggregation-unit-test-v1",
    "status": "pass_aggregation_structure_unit_tests",
    "tests_passed": 5,
    "control_prediction_member_opened": False,
    "outcomes_read": False,
    "metrics_calculated": False,
}
(root / "test_receipt.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
freeze_tree(root, {"artifact_class": "lsgkm_aggregation_structure_tests", **receipt})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Aggregation tests: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
