#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-044
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-044-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-044-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
EVALUATOR=${ROOT}/src/masld_bench/evaluators/rna_atac_cobolt_5seed_v8_development.py
TEST=${ROOT}/tests/unit/test_rna_atac_cobolt_5seed_v8_evaluator.py
VERIFIER=${ROOT}/slurm/verify_cobolt_v8_cpu_baselines.sbatch
SOURCE_JOB=${ROOT}/slurm/test_and_freeze_rna_atac_cobolt_5seed_v8_evaluator.sbatch
BUILDER=${ROOT}/executions/dispatch-builders/cobolt-5seed-v8/build_bundle.py
BUILD_JOB=${ROOT}/executions/dispatch-builders/cobolt-5seed-v8/build_bundle.sbatch
CANARY=${ROOT}/executions/dispatch-bundles/cobolt-5seed-cpu-submit-v1/submit_cpu_baseline_canary_v8.sbatch
REMAINING=${ROOT}/executions/dispatch-bundles/cobolt-5seed-cpu-submit-v1/submit_cpu_baselines_remaining_v8.sbatch
EVALUATION=${ROOT}/slurm/evaluate_rna_atac_cobolt_5seed_v8_smoke.sbatch
TARGET=${ROOT}/executions/cobolt-v8-tooling-test-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite v8 tooling test: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
bash -n "${VERIFIER}" "${SOURCE_JOB}" "${BUILD_JOB}" "${CANARY}" "${REMAINING}" "${EVALUATION}"
PYTHONPATH=src python -m py_compile "${EVALUATOR}" "${TEST}" "${BUILDER}"
PYTHONPATH=src python -m unittest \
  tests.unit.test_rna_atac_cobolt_5seed_v8_evaluator \
  tests.unit.test_rna_atac_cobolt_5seed_evaluator \
  tests.unit.test_rna_atac_cobolt_evaluator
PYTHONPATH=src python - "${EVALUATOR}" "${VERIFIER}" "${SOURCE_JOB}" "${BUILDER}" "${BUILD_JOB}" "${CANARY}" "${REMAINING}" "${EVALUATION}" "${STAGE}" <<'PY'
from pathlib import Path
import importlib.util
import json
import sys

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.hashing import sha256_file

evaluator, verifier, source_job, builder, build_job, canary, remaining, evaluation, stage = map(Path, sys.argv[1:])
spec = importlib.util.spec_from_file_location("cobolt_v8_builder", builder)
if spec is None or spec.loader is None:
    raise SystemExit("could not import v8 builder")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
if module.CAMPAIGN_ID != "v1-rna-atac-cobolt-5seed-smoke-v8":
    raise SystemExit("builder campaign binding differs")
if module.BUNDLE_ID != "cobolt-5seed-25run-v8-generic5":
    raise SystemExit("builder bundle binding differs")
if module.WRAPPER_NAME != "model_work_205.sbatch":
    raise SystemExit("builder wrapper name differs")
for path in (verifier, source_job, build_job, canary, remaining, evaluation):
    text = path.read_text(encoding="utf-8")
    for required in (
        "#SBATCH --partition=cpu",
        "#SBATCH --qos=nslab",
        "COBOLT_V8_CANDIDATE",
        "COBOLT_V8_CANDIDATE_SHA",
        "v1-rna-atac-cobolt-5seed-smoke-v8--",
    ):
        if required not in text:
            raise SystemExit(f"{path.name} lacks {required}")
builder_text = builder.read_text(encoding="utf-8")
for required in (
    "#SBATCH --job-name=model-work-205",
    "#SBATCH --account=nslab",
    "#SBATCH --partition=gpu",
    "#SBATCH --qos=nslab",
    "#SBATCH --gres=gpu:l40s:1",
    '"cpus": 4',
    '"memory_gb": 64',
    '"walltime": "06:00:00"',
):
    if required not in builder_text:
        raise SystemExit(f"builder lacks {required}")
records = {
    path.name: {"sha256": sha256_file(path), "size_bytes": path.stat().st_size}
    for path in (evaluator, verifier, source_job, builder, build_job, canary, remaining, evaluation)
}
write_json_exclusive(stage / "tooling_test_receipt.json", {
    "schema_version": "masld-bench-cobolt-v8-tooling-test-v1",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v8",
    "records": records,
    "unit_tests_passed": 9,
    "tooling_files_statically_checked": 8,
    "outcomes_opened": False,
    "gpu_submitted": False,
    "status": "pass",
})
print(freeze_tree(stage, {
    "artifact_class": "cobolt_v8_tooling_test",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v8",
    "outcomes_opened": False,
    "gpu_submitted": False,
    "status": "pass",
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
