#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-verify-207
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-verify-207-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-verify-207-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v7--a717faf7669b6f26
CAMPAIGN_SHA=2dedff948bdb8f97cda296ddf8b56329e738185ab17939bf2fad2682cf1bf85d
EXECUTION_ROOT=${ROOT}/candidates/executions/v1-rna-atac-cobolt-5seed-smoke-v7--a717faf7669b6f26
TARGET=${ROOT}/executions/cobolt-v7-cpu-baseline-verification-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite CPU verification: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
printf '%s  %s\n' "${CAMPAIGN_SHA}" "${CANDIDATE}/ARTIFACTS.json" \
  | sha256sum --check --strict
PYTHONPATH=src python - \
  "${CANDIDATE}" "${EXECUTION_ROOT}" "${STAGE}" "${CAMPAIGN_SHA}" <<'PY'
from pathlib import Path
import csv
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.campaign import verify_run_execution_attempt
from masld_bench.hashing import sha256_file
from masld_bench.planner import load_frozen_plan

candidate = Path(sys.argv[1]).resolve(strict=True)
execution = Path(sys.argv[2]).resolve(strict=True)
stage = Path(sys.argv[3]).resolve(strict=True)
candidate_sha256 = sys.argv[4]
verify_frozen_tree(candidate)
plan = load_frozen_plan(candidate)
runs = sorted(
    (
        run
        for run in plan["runs"]
        if run.get("resource_profile") == "cpu_baseline_smoke"
    ),
    key=lambda run: str(run["run_id"]),
)
if len(runs) != 150:
    raise SystemExit(f"expected 150 CPU baseline runs, observed {len(runs)}")
records = []
for index, run in enumerate(runs, start=1):
    attempt = execution / "runs" / str(run["run_id"]) / "attempt-001"
    receipt = verify_run_execution_attempt(attempt, require_succeeded=True)
    records.append({
        "index": index,
        "run_id": str(run["run_id"]),
        "model_id": str(run["model_id"]),
        "seed": int(run["seed"]),
        "fold": int(run["fold"]),
        "status": str(receipt["status"]),
        "attempt_path": attempt.as_posix(),
        "attempt_artifacts_sha256": sha256_file(attempt / "ARTIFACTS.json"),
    })
fields = (
    "index",
    "run_id",
    "model_id",
    "seed",
    "fold",
    "status",
    "attempt_path",
    "attempt_artifacts_sha256",
)
with (stage / "verified_attempts.tsv").open("x", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
write_json_exclusive(stage / "verification_receipt.json", {
    "schema_version": "masld-bench-cobolt-v7-cpu-verification-v1",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v7",
    "campaign_artifact_sha256": candidate_sha256,
    "expected_attempts": 150,
    "verified_succeeded_attempts": len(records),
    "outcome_records_parsed_or_joined": False,
    "verification_scope": "immutable run, adapter-output, and prediction artifacts",
    "status": "pass",
})
artifact_hash = freeze_tree(stage, {
    "artifact_class": "cobolt_v7_cpu_baseline_verification",
    "campaign_sha256": candidate_sha256,
    "verified_attempts": len(records),
    "status": "pass",
})
print(artifact_hash)
PY
mv "${STAGE}" "${TARGET}"
printf 'Cobolt v7 CPU verification: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
