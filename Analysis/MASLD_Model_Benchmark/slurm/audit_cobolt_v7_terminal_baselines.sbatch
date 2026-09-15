#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-audit-207
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-207-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-207-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v7--a717faf7669b6f26
CANDIDATE_SHA=2dedff948bdb8f97cda296ddf8b56329e738185ab17939bf2fad2682cf1bf85d
EXECUTION_ROOT=${ROOT}/candidates/executions/v1-rna-atac-cobolt-5seed-smoke-v7--a717faf7669b6f26
TARGET=${ROOT}/executions/cobolt-v7-terminal-baseline-audit-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite terminal audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
printf '%s  %s\n' "${CANDIDATE_SHA}" "${CANDIDATE}/ARTIFACTS.json" | sha256sum --check --strict
PYTHONPATH=src python - "${CANDIDATE}" "${EXECUTION_ROOT}" "${STAGE}" "${CANDIDATE_SHA}" <<'PY'
from pathlib import Path
from collections import Counter
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file
from masld_bench.planner import load_frozen_plan

candidate = Path(sys.argv[1]).resolve(strict=True)
execution = Path(sys.argv[2]).resolve(strict=True)
stage = Path(sys.argv[3]).resolve(strict=True)
candidate_sha256 = sys.argv[4]
verify_frozen_tree(candidate)
plan = load_frozen_plan(candidate)
runs = sorted(
    (run for run in plan["runs"] if run.get("resource_profile") == "cpu_baseline_smoke"),
    key=lambda run: str(run["run_id"]),
)
if len(runs) != 150:
    raise SystemExit(f"expected 150 CPU baselines, observed {len(runs)}")
records = []
for run in runs:
    attempt = execution / "runs" / str(run["run_id"]) / "attempt-001"
    verify_frozen_tree(attempt)
    receipt = json.loads((attempt / "run_execution_receipt.json").read_text(encoding="utf-8"))
    status = receipt.get("status")
    actions = receipt.get("action_receipts")
    if (
        receipt.get("campaign_id") != "v1-rna-atac-cobolt-5seed-smoke-v7"
        or receipt.get("candidate_manifest_sha256") != candidate_sha256
        or receipt.get("model_id") != run["model_id"]
        or receipt.get("seed") != run["seed"]
        or receipt.get("fold") != run["fold"]
        or status not in {"succeeded", "invalid_output"}
        or not isinstance(actions, list)
        or [item.get("action") for item in actions] != ["prepare", "fit", "predict"]
        or any(item.get("status") != "complete" for item in actions)
    ):
        raise SystemExit(f"terminal attempt contract differs: {attempt}")
    failure = receipt.get("failure")
    if status == "succeeded":
        if failure is not None or receipt.get("retry_class") is not None:
            raise SystemExit(f"successful terminal receipt differs: {attempt}")
    elif (
        receipt.get("retry_class") != "software"
        or not isinstance(failure, str)
        or "FirewallError: current MASLD Resource authorities changed during model work" not in failure
        or "docs/" not in failure
    ):
        raise SystemExit(f"invalid terminal failure differs: {attempt}")
    records.append({
        "run_id": str(run["run_id"]),
        "model_id": str(run["model_id"]),
        "seed": int(run["seed"]),
        "fold": int(run["fold"]),
        "status": str(status),
        "retry_class": "" if receipt.get("retry_class") is None else str(receipt["retry_class"]),
        "attempt_artifacts_sha256": sha256_file(attempt / "ARTIFACTS.json"),
        "adapter_actions_complete": True,
        "failure": "" if failure is None else failure,
    })
counts = Counter(record["status"] for record in records)
if counts != {"succeeded": 42, "invalid_output": 108}:
    raise SystemExit(f"terminal disposition census differs: {dict(counts)}")
fields = (
    "run_id", "model_id", "seed", "fold", "status", "retry_class",
    "attempt_artifacts_sha256", "adapter_actions_complete", "failure",
)
with (stage / "terminal_attempts.tsv").open("x", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
write_json_exclusive(stage / "audit_receipt.json", {
    "schema_version": "masld-bench-cobolt-v7-terminal-baseline-audit-v1",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v7",
    "candidate_artifacts_sha256": candidate_sha256,
    "expected_attempts": 150,
    "verified_attempts": len(records),
    "succeeded": counts["succeeded"],
    "invalid_output_resource_firewall": counts["invalid_output"],
    "all_adapter_actions_complete": True,
    "outcomes_parsed_or_joined": False,
    "gpu_gate_passed": False,
    "status": "campaign_invalidated",
})
print(freeze_tree(stage, {
    "artifact_class": "cobolt_v7_terminal_baseline_audit",
    "candidate_artifacts_sha256": candidate_sha256,
    "verified_attempts": len(records),
    "succeeded": counts["succeeded"],
    "invalid_output": counts["invalid_output"],
    "gpu_gate_passed": False,
    "status": "campaign_invalidated",
}))
PY
mv "${STAGE}" "${TARGET}"
printf 'Cobolt v7 terminal baseline audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
