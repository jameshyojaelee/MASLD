#!/usr/bin/env bash
#SBATCH --job-name=masld-seq-bundle001
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=20:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${MASLD_GPU_DISPATCHER_AUTHORITY:?dispatcher authority is required}"
if [[ "${MASLD_GPU_DISPATCHER_AUTHORITY}" != sequence_gpu_dispatcher_v2 ]]; then
  printf 'Invalid GPU dispatcher authority.\n' >&2
  exit 2
fi
: "${SLURM_JOB_ID:?SLURM_JOB_ID is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
MANIFEST="${ROOT}/config/campaigns/sequence_gpu_bundle_001_20260824.tsv"
MANIFEST_SHA256=042faacf41d698e4003e5ed18829530aa7e9c69b8d3cc5c375775027b7bfcad1
TARGET="${ROOT}/executions/sequence-gpu-bundle-001-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

printf '%s  %s\n' "${MANIFEST_SHA256}" "${MANIFEST}" | sha256sum --check --strict
if [[ -d "${TARGET}" ]]; then
  printf 'Bundle allocation already terminal: %s\n' "${TARGET}"
  exit 0
fi
install -d -m 0750 "${STAGE}/logs" "${STAGE}/task_receipts"

while IFS=$'\t' read -r task_order logical_work_id model_id stage lineage_id \
  outer_fold seed input_job_id input_artifacts_sha256 script \
  expected_artifact_template prerequisite_logical_id failure_scope; do
  [[ "${task_order}" == task_order ]] && continue
  receipt="${STAGE}/task_receipts/${logical_work_id}.tsv"
  if [[ -f "${receipt}" ]]; then
    printf 'Task checkpoint already exists, skipping: %s\n' "${logical_work_id}"
    continue
  fi
  if [[ "${prerequisite_logical_id}" != none ]]; then
    printf 'Unexpected prerequisite in bundle 001: %s\n' "${logical_work_id}" >&2
    exit 2
  fi
  if [[ "${failure_scope}" != isolated_continue ]]; then
    printf 'Unexpected failure scope in bundle 001: %s\n' "${logical_work_id}" >&2
    exit 2
  fi
  expected_relative="${expected_artifact_template//\{SLURM_JOB_ID\}/${SLURM_JOB_ID}}"
  expected="${ROOT}/${expected_relative}"
  started_at="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  status=failed
  return_code=1
  artifact_sha256=not_available
  if [[ -f "${expected}/ARTIFACTS.json" ]]; then
    status=recovered_passed
    return_code=0
    read -r artifact_sha256 _ < <(sha256sum "${expected}/ARTIFACTS.json")
  elif [[ -d "${expected}.failed" ]]; then
    status=existing_failed
  else
    set +e
    env \
      OUTER_FOLD="${outer_fold}" \
      LINEAGE="${lineage_id}" \
      SEED="${seed}" \
      INPUT_JOB_ID="${input_job_id}" \
      INPUT_ARTIFACTS_SHA256="${input_artifacts_sha256}" \
      bash "${ROOT}/${script}" \
      > "${STAGE}/logs/${logical_work_id}.stdout.txt" \
      2> "${STAGE}/logs/${logical_work_id}.stderr.txt"
    return_code=$?
    set -e
    if (( return_code == 0 )) && [[ -f "${expected}/ARTIFACTS.json" ]]; then
      status=passed
      read -r artifact_sha256 _ < <(sha256sum "${expected}/ARTIFACTS.json")
    elif (( return_code == 0 )); then
      status=failed_missing_artifact
      return_code=97
    fi
  fi
  ended_at="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  receipt_tmp="${receipt}.tmp"
  printf 'logical_work_id\tmodel_id\tlineage_id\touter_fold\tseed\tstatus\treturn_code\texpected_artifact\tartifacts_sha256\tstarted_at\tended_at\n' \
    > "${receipt_tmp}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${logical_work_id}" "${model_id}" "${lineage_id}" "${outer_fold}" \
    "${seed}" "${status}" "${return_code}" "${expected}" \
    "${artifact_sha256}" "${started_at}" "${ended_at}" >> "${receipt_tmp}"
  mv "${receipt_tmp}" "${receipt}"
  printf 'Task %s terminal status: %s\n' "${logical_work_id}" "${status}"
done < "${MANIFEST}"

sha256sum \
  "${MANIFEST}" \
  "${ROOT}/slurm/sequence_gpu_bundle_001.sbatch" \
  "${ROOT}/slurm/bpnet_train_predict_hepatocyte_donor0_genomic0.sbatch" \
  "${ROOT}/slurm/train_chrombpnet_hepatocyte_donor0_genomic0_seed20260824.sbatch" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${MANIFEST}" <<'PY'
from pathlib import Path
from collections import Counter
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
manifest_path = Path(sys.argv[2]).resolve(strict=True)
with manifest_path.open(encoding="utf-8", newline="") as handle:
    tasks = list(csv.DictReader(handle, delimiter="\t"))
receipts = []
for task in tasks:
    path = root / "task_receipts" / f"{task['logical_work_id']}.tsv"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1 or rows[0]["logical_work_id"] != task["logical_work_id"]:
        raise SystemExit("bundle task receipt differs")
    receipt = rows[0]
    if receipt["status"] in {"passed", "recovered_passed"}:
        artifact = Path(receipt["expected_artifact"])
        verify_frozen_tree(artifact)
        actual = __import__("hashlib").sha256(
            (artifact / "ARTIFACTS.json").read_bytes()
        ).hexdigest()
        if actual != receipt["artifacts_sha256"]:
            raise SystemExit("task artifact receipt hash differs")
    receipts.append(receipt)
counts = Counter(receipt["status"] for receipt in receipts)
terminal_status = (
    "passed" if set(counts) <= {"passed", "recovered_passed"}
    else "completed_with_task_failures"
)
(root / "bundle_summary.json").write_text(
    json.dumps({
        "schema_version": "masld-bench-sequence-gpu-bundle-receipt-v1",
        "bundle_id": "sequence_gpu_bundle_001",
        "logical_tasks": len(tasks),
        "status_counts": dict(sorted(counts.items())),
        "terminal_status": terminal_status,
        "continue_after_isolated_failure": True,
        "checkpoint_granularity": "logical_task",
    }, sort_keys=True, indent=2) + "\n",
    encoding="utf-8",
)
freeze_tree(root, {
    "artifact_class": "sequence_gpu_bundle_allocation_receipt",
    "bundle_id": "sequence_gpu_bundle_001",
    "logical_tasks": len(tasks),
    "terminal_status": terminal_status,
    "continue_after_isolated_failure": True,
    "checkpoint_granularity": "logical_task",
})
PY
module purge
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
