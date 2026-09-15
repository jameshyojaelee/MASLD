#!/usr/bin/env bash
#SBATCH --job-name=model-check-214
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-check-214-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-check-214-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${COBOLT_V9_CANDIDATE:?set COBOLT_V9_CANDIDATE to the exact frozen production candidate}
CANDIDATE_SHA=${COBOLT_V9_CANDIDATE_SHA:?set COBOLT_V9_CANDIDATE_SHA to its ARTIFACTS.json sha256}
STABILITY=${COBOLT_V9_STABILITY:?set COBOLT_V9_STABILITY to the exact frozen stability artifact}
STABILITY_SHA=${COBOLT_V9_STABILITY_SHA:?set COBOLT_V9_STABILITY_SHA to its ARTIFACTS.json sha256}
EXECUTOR=${ROOT}/scripts/run_cobolt_v9_cpu_baseline_bundle.py
EXECUTOR_SHA=${COBOLT_V9_CPU_EXECUTOR_SHA:?set COBOLT_V9_CPU_EXECUTOR_SHA to the reviewed executor sha256}
TARGET=${ROOT}/executions/cobolt-v9-cpu-bundle-validation-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite validation output: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/model-check-214-${SLURM_JOB_ID}
cd "${ROOT}"
printf '%s  %s\n' "${CANDIDATE_SHA}" "${CANDIDATE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${STABILITY_SHA}" "${STABILITY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${EXECUTOR_SHA}" "${EXECUTOR}" | sha256sum --check --strict
PYTHONPATH=src python -m py_compile "${EXECUTOR}"
PYTHONPATH=src python -m py_compile tests/unit/test_cobolt_v9_cpu_baseline_bundle.py
bash -n slurm/run_cobolt_v9_current_cpu_baseline_bundle.sbatch
PYTHONPATH=src python -m unittest tests.unit.test_cobolt_v9_cpu_baseline_bundle \
  >"${STAGE}/tests.stdout.txt" 2>"${STAGE}/tests.stderr.txt"
PYTHONPATH=src python "${EXECUTOR}" \
  --candidate "${CANDIDATE}" \
  --candidate-sha256 "${CANDIDATE_SHA}" \
  --stability "${STABILITY}" \
  --stability-sha256 "${STABILITY_SHA}" \
  --executor-sha256 "${EXECUTOR_SHA}" \
  --workers 12 \
  --allocated-cpus 16 \
  --allocated-memory-gb 200 \
  --preflight >"${STAGE}/preflight.json"
sha256sum \
  "${EXECUTOR}" \
  tests/unit/test_cobolt_v9_cpu_baseline_bundle.py \
  slurm/run_cobolt_v9_current_cpu_baseline_bundle.sbatch \
  slurm/validate_cobolt_v9_cpu_baseline_bundle.sbatch \
  >"${STAGE}/source.sha256"

PYTHONPATH=src python - "${STAGE}" "${CANDIDATE_SHA}" "${STABILITY_SHA}" "${EXECUTOR_SHA}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, write_json_exclusive

stage = Path(sys.argv[1]).resolve(strict=True)
preflight = json.loads((stage / "preflight.json").read_text(encoding="utf-8"))
if (
    preflight.get("status") != "pass"
    or preflight.get("runs") != 150
    or preflight.get("execution_started") is not False
    or preflight.get("sealed_features_opened") is not False
    or preflight.get("sealed_labels_opened") is not False
    or preflight.get("outcome_records_parsed_or_joined") is not False
):
    raise SystemExit("CPU bundle preflight differs")
write_json_exclusive(stage / "validation_receipt.json", {
    "schema_version": "masld-bench-cobolt-v9-cpu-bundle-validation-v1",
    "candidate_manifest_sha256": sys.argv[2],
    "stability_manifest_sha256": sys.argv[3],
    "executor_sha256": sys.argv[4],
    "focused_tests": 7,
    "preflight_runs": 150,
    "execution_started": False,
    "outcomes_opened": False,
    "status": "pass",
})
print(freeze_tree(stage, {
    "artifact_class": "cobolt_v9_cpu_bundle_validation",
    "candidate_manifest_sha256": sys.argv[2],
    "stability_manifest_sha256": sys.argv[3],
    "executor_sha256": sys.argv[4],
    "status": "pass",
    "outcomes_opened": False,
}))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
