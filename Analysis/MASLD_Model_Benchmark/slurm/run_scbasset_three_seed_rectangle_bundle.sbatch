#!/usr/bin/env bash
#SBATCH --job-name=model-training
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

: "${SCBASSET_RECTANGLE_AUTHORITY:?central dispatcher authority is required}"
: "${SCBASSET_RECTANGLE_CONTRACT:?rectangle contract is required}"
: "${SCBASSET_RECTANGLE_CONTRACT_SHA256:?rectangle contract hash is required}"
: "${SCBASSET_PREFLIGHT_ONLY:=0}"
if [[ "${SCBASSET_RECTANGLE_AUTHORITY}" != scbasset_three_seed_rectangle_v1 ]]; then
  printf 'Invalid scBasset rectangle authority\n' >&2
  exit 2
fi
if [[ "${SCBASSET_PREFLIGHT_ONLY}" != 0 && "${SCBASSET_PREFLIGHT_ONLY}" != 1 ]]; then
  printf 'Invalid scBasset preflight mode\n' >&2
  exit 2
fi

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
EXPECTED_CONTRACT="${ROOT}/config/campaigns/scbasset_three_seed_rectangle_20260825.json"
TRAIN="${ROOT}/slurm/scbasset_train_donor0_genomic0_seed11.sbatch"
PREDICT="${ROOT}/slurm/scbasset_predict_valid_donor0_genomic0_seed11.sbatch"
FAILURE_POLICY=continue_after_isolated_fit_failure
if [[ "${SCBASSET_RECTANGLE_CONTRACT}" != "${EXPECTED_CONTRACT}" ]]; then
  printf 'scBasset rectangle contract path differs\n' >&2
  exit 2
fi
printf '%s  %s\n' "${SCBASSET_RECTANGLE_CONTRACT_SHA256}" "${EXPECTED_CONTRACT}" \
  | sha256sum --check --strict

if [[ "${SCBASSET_PREFLIGHT_ONLY}" == 1 ]]; then
  : "${SCBASSET_PREFLIGHT_STAGE:?CPU preflight stage is required}"
  case "${SCBASSET_PREFLIGHT_STAGE}" in
    "${ROOT}"/executions/scbasset-three-seed-rectangle-admission-20260825.staging-*/validation/bundle-preflight) ;;
    *) printf 'Preflight stage is outside the admitted namespace\n' >&2; exit 2 ;;
  esac
  STAGE="${SCBASSET_PREFLIGHT_STAGE}"
  install -d -m 0750 "${STAGE}"
else
  : "${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
  TARGET="${ROOT}/executions/scbasset-three-seed-rectangle-${SLURM_JOB_ID}"
  STAGE="${TARGET}.staging"
  FAILED="${TARGET}.failed"
  preserve_failed_stage() {
    if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
      mv "${STAGE}" "${FAILED}"
    fi
  }
  trap preserve_failed_stage EXIT
  if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
    printf 'Refusing to overwrite scBasset rectangle bundle: %s\n' "${TARGET}" >&2
    exit 1
  fi
  install -d -m 0750 "${STAGE}/logs" "${STAGE}/fit_receipts"
fi

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${ROOT}" "${EXPECTED_CONTRACT}" "${STAGE}" \
  "${SCBASSET_PREFLIGHT_ONLY}" <<'PY'
from hashlib import sha256
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
contract_path = Path(sys.argv[2]).resolve(strict=True)
stage = Path(sys.argv[3]).resolve(strict=True)
preflight_only = sys.argv[4] == "1"
contract = json.loads(contract_path.read_text(encoding="utf-8"))
if (
    contract.get("schema_version") != "masld-bench-scbasset-three-seed-rectangle-v1"
    or contract.get("status")
    != "PROSPECTIVE_DEVELOPMENT_RECTANGLE_WAITING_FOR_CENTRAL_DISPATCH"
    or contract.get("production_bundle", {}).get("bundle_id") != "model-training-604"
    or contract.get("production_bundle", {}).get("fits") != 11
    or contract.get("production_bundle", {}).get("failure_policy")
    != "continue_after_isolated_fit_failure"
):
    raise SystemExit("scBasset rectangle contract differs")
sources = contract["source_authorities"]
for record in sources.values():
    path = root / record["path"]
    if sha256(path.read_bytes()).hexdigest() != record["sha256"]:
        raise SystemExit(f"source authority differs: {path}")
admission_root = root / contract["production_bundle"]["admission_artifact_root"]
if not preflight_only:
    admission = verify_frozen_tree(admission_root)
    if (
        admission.get("metadata", {}).get("artifact_class")
        != "scbasset_three_seed_rectangle_admission"
        or admission.get("metadata", {}).get("contract_sha256")
        != sha256(contract_path.read_bytes()).hexdigest()
        or admission.get("metadata", {}).get("status") != "passed"
    ):
        raise SystemExit("scBasset admission gate differs")
inputs = {int(record["outer_fold"]): record for record in contract["input_artifacts"]}
missing = contract["missing_contexts"]
if len(inputs) != 5 or len(missing) != 11:
    raise SystemExit("scBasset input or missing-fit census differs")
with (stage / "fit_plan.tsv").open("x", encoding="utf-8", newline="") as handle:
    fields = (
        "fit_order",
        "outer_fold",
        "seed",
        "split_id",
        "input_job_id",
        "input_root",
        "input_artifacts_sha256",
    )
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for order, (fold, seed) in enumerate(missing, start=1):
        record = inputs[int(fold)]
        writer.writerow(
            {
                "fit_order": order,
                "outer_fold": fold,
                "seed": seed,
                "split_id": f"donor{fold}_genomic{fold}",
                "input_job_id": record["job_id"],
                "input_root": str(root / record["artifact_root"]),
                "input_artifacts_sha256": record["artifacts_sha256"],
            }
        )
(stage / "wrapper_preflight.json").write_text(
    json.dumps(
        {
            "status": "pass",
            "bundle_id": "model-training-604",
            "fits": 11,
            "logical_tasks": 22,
            "failure_policy": "continue_after_isolated_fit_failure",
            "outcome_paths_available": False,
            "sealed_paths_available": False,
            "direct_gpu_submission_performed": False,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
PY
bash -n "${TRAIN}"
bash -n "${PREDICT}"
if [[ "${SCBASSET_PREFLIGHT_ONLY}" == 1 ]]; then
  module purge
  printf 'scBasset rectangle CPU preflight passed\n'
  exit 0
fi
module purge

printf 'fit_order\touter_fold\tseed\ttrain_rc\tpredict_rc\tmodel_path\tprediction_path\n' \
  > "${STAGE}/dispositions.tsv"
failed=0
while IFS=$'\t' read -r fit_order outer_fold seed split_id input_job_id input_root input_sha; do
  [[ "${fit_order}" == fit_order ]] && continue
  model_path="${ROOT}/executions/scbasset-${split_id}-seed${seed}-${SLURM_JOB_ID}"
  prediction_path="${ROOT}/executions/scbasset-${split_id}-seed${seed}-valid-${SLURM_JOB_ID}"
  train_rc=0
  predict_rc=NA
  if OUTER_FOLD="${outer_fold}" SEED="${seed}" INPUT_JOB_ID="${input_job_id}" \
      INPUT_ARTIFACTS_SHA256="${input_sha}" INPUT_ROOT_OVERRIDE="${input_root}" \
      bash "${TRAIN}" > "${STAGE}/logs/fit${fit_order}.train.stdout.txt" \
      2> "${STAGE}/logs/fit${fit_order}.train.stderr.txt"; then
    if OUTER_FOLD="${outer_fold}" SEED="${seed}" INPUT_JOB_ID="${input_job_id}" \
        INPUT_ARTIFACTS_SHA256="${input_sha}" INPUT_ROOT_OVERRIDE="${input_root}" \
        MODEL_JOB_ID="${SLURM_JOB_ID}" MODEL_ARTIFACTS_SHA256=AUTO RUN_ID=AUTO \
        bash "${PREDICT}" > "${STAGE}/logs/fit${fit_order}.predict.stdout.txt" \
        2> "${STAGE}/logs/fit${fit_order}.predict.stderr.txt"; then
      predict_rc=0
    else
      predict_rc=$?
      failed=1
    fi
  else
    train_rc=$?
    failed=1
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${fit_order}" "${outer_fold}" "${seed}" "${train_rc}" "${predict_rc}" \
    "${model_path}" "${prediction_path}" >> "${STAGE}/dispositions.tsv"
done < "${STAGE}/fit_plan.tsv"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${SCBASSET_RECTANGLE_CONTRACT_SHA256}" \
  "${SLURM_JOB_ID}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.hashing import sha256_file

stage = Path(sys.argv[1]).resolve(strict=True)
contract_sha256 = sys.argv[2]
job_id = sys.argv[3]
rows = list(csv.DictReader((stage / "dispositions.tsv").open(), delimiter="\t"))
receipts = []
for row in rows:
    model = Path(row["model_path"])
    prediction = Path(row["prediction_path"])
    model_hash = "not_available"
    prediction_hash = "not_available"
    if row["train_rc"] == "0" and (model / "ARTIFACTS.json").is_file():
        verify_frozen_tree(model)
        model_hash = sha256_file(model / "ARTIFACTS.json")
    if row["predict_rc"] == "0" and (prediction / "ARTIFACTS.json").is_file():
        verify_frozen_tree(prediction)
        prediction_hash = sha256_file(prediction / "ARTIFACTS.json")
    receipt = {
        "fit_order": int(row["fit_order"]),
        "outer_fold": int(row["outer_fold"]),
        "seed": int(row["seed"]),
        "train_returncode": int(row["train_rc"]),
        "predict_returncode": None if row["predict_rc"] == "NA" else int(row["predict_rc"]),
        "model_path": str(model),
        "model_artifacts_sha256": model_hash,
        "prediction_path": str(prediction),
        "prediction_artifacts_sha256": prediction_hash,
    }
    (stage / "fit_receipts" / f"fit{int(row['fit_order']):02d}.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipts.append(receipt)
status = "succeeded" if all(
    row["train_returncode"] == 0 and row["predict_returncode"] == 0
    for row in receipts
) else "partial_failure"
(stage / "bundle_summary.json").write_text(
    json.dumps(
        {
            "schema_version": "masld-bench-scbasset-three-seed-bundle-receipt-v1",
            "status": status,
            "bundle_id": "model-training-604",
            "slurm_job_id": job_id,
            "fits": 11,
            "logical_tasks": 22,
            "dispositions": receipts,
            "failure_policy": "continue_after_isolated_fit_failure",
            "prediction_values_read_by_controller": False,
            "metric_values_read_by_controller": False,
            "raw_outcomes_read": False,
            "sealed_data_read": False,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
freeze_tree(
    stage,
    {
        "artifact_class": "scbasset_three_seed_rectangle_bundle",
        "bundle_id": "model-training-604",
        "slurm_job_id": job_id,
        "fits": 11,
        "logical_tasks": 22,
        "terminal_status": status,
        "failure_policy": "continue_after_isolated_fit_failure",
        "contract_sha256": contract_sha256,
        "prediction_values_read": False,
        "metric_values_read": False,
        "raw_outcomes_read": False,
        "sealed_data_read": False,
        "status": status,
    },
)
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'scBasset three-seed rectangle bundle: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
exit "${failed}"
