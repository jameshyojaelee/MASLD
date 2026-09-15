#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-401
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
plan="${root}/executions/sequence-gpu-bundle-plan-21068430"
wrapper="${root}/slurm/run_sequence_gpu_bundle.sbatch"
target="${root}/executions/sequence-gpu-bundles-next-validation-${SLURM_JOB_ID}"
stage="${target}.staging"
failed="${target}.failed"
plan_sha256="36dd9ad437bdd61e150e1c681a88bfdf946672c768bffe5f78fac8ab2c352e69"
wrapper_sha256="b1f5bd6113d937c009feaccd85b4351e3bc05aeaee9b2f17ac12c980168c162d"
bundle_ids=(
  sequence_gpu_bundle_016_sequence_cnn_control-ready_36823fe74c
  sequence_gpu_bundle_023_sequence_transformer_control-ready_ea978e9249
  sequence_gpu_bundle_004_bpnet-ready_71d22f6453
  sequence_gpu_bundle_010_chrombpnet-ready_a7f761b0ef
)

preserve_failed_stage() {
  if [[ -d "${stage}" && ! -e "${failed}" ]]; then
    mv "${stage}" "${failed}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${target}" || -e "${stage}" || -e "${failed}" ]]; then
  printf 'Refusing to overwrite next-bundle validation artifact.\n' >&2
  exit 1
fi
install -d -m 0750 "${stage}/validation"
printf '%s  %s\n' "${plan_sha256}" "${plan}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${wrapper_sha256}" "${wrapper}" | sha256sum --check --strict
bash -n "${wrapper}"

for bundle_id in "${bundle_ids[@]}"; do
  install -d -m 0750 \
    "${stage}/validation/${bundle_id}/logs" \
    "${stage}/validation/${bundle_id}/task_receipts"
  env \
    MASLD_GPU_DISPATCHER_AUTHORITY=sequence_gpu_dispatcher_v2 \
    BUNDLE_PLAN_ROOT="${plan}" \
    BUNDLE_PLAN_ARTIFACTS_SHA256="${plan_sha256}" \
    BUNDLE_ID="${bundle_id}" \
    BUNDLE_REQUESTED_TIME=40:00:00 \
    BUNDLE_PREFLIGHT_ONLY=1 \
    BUNDLE_PREFLIGHT_STAGE="${stage}/validation/${bundle_id}" \
    bash "${wrapper}"
done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${root}/src" python - "${stage}" "${plan}" "${bundle_ids[@]}" <<'PY'
from pathlib import Path
from collections import Counter
import csv
import json
import sys
from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
plan = Path(sys.argv[2]).resolve(strict=True)
bundle_ids = sys.argv[3:]
with (plan / "plan/bundle_allocations.tsv").open(encoding="utf-8", newline="") as handle:
    allocations = {row["bundle_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
with (plan / "plan/bundle_tasks.tsv").open(encoding="utf-8", newline="") as handle:
    tasks = list(csv.DictReader(handle, delimiter="\t"))
expected_models = {
    bundle_ids[0]: "sequence_cnn_control",
    bundle_ids[1]: "sequence_transformer_control",
    bundle_ids[2]: "bpnet",
    bundle_ids[3]: "chrombpnet",
}
seen = set()
for bundle_id in bundle_ids:
    allocation = allocations[bundle_id]
    selected = [row for row in tasks if row["bundle_id"] == bundle_id]
    preflight = json.loads(
        (stage / "validation" / bundle_id / "preflight.json").read_text()
    )
    if (
        len(selected) != 8
        or allocation["initial_readiness"] != "ready"
        or allocation["requested_time"] != "40:00:00"
        or allocation["partition"] != "gpu"
        or allocation["qos"] != "nslab"
        or Counter(row["model_id"] for row in selected) != {expected_models[bundle_id]: 8}
        or any(row["current_state"] != "backlog_ready" for row in selected)
        or preflight["status"] != "pass"
        or preflight["logical_tasks"] != 8
        or preflight["outcome_paths_available_to_models"]
        or preflight["sealed_paths_available_to_models"]
    ):
        raise SystemExit(f"next sequence bundle contract differs: {bundle_id}")
    logical_ids = {row["logical_work_id"] for row in selected}
    if seen & logical_ids:
        raise SystemExit("next sequence bundles overlap")
    seen |= logical_ids
summary = {
    "schema_version": "masld-bench-sequence-gpu-bundle-validation-v3",
    "status": "pass",
    "bundles": bundle_ids,
    "logical_tasks": len(seen),
    "gpu_jobs_submitted": 0,
    "future_qos": "nslab",
    "generic_job_name": "model-training-400",
    "sealed_paths_available_to_models": False,
}
(stage / "validation/summary.json").write_text(
    json.dumps(summary, sort_keys=True, indent=2) + "\n"
)
freeze_tree(stage, {
    "artifact_class": "sequence_gpu_bundle_validation",
    **summary,
})
PY

trap - EXIT
mv "${stage}" "${target}"
sha256sum "${target}/ARTIFACTS.json"
