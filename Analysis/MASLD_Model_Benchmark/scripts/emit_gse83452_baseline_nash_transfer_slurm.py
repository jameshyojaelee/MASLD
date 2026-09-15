#!/usr/bin/env python3
"""Emit the downstream GSE83452 transfer sbatch scripts from frozen upstream hashes.

The registration, prediction, evaluator-label and evaluation steps each pin the
SHA-256 of the output-file tree they consume, and those hashes only exist once the
upstream job has landed.  Writing the sbatch files by hand at that point is how a
hash gets transcribed wrong, so they are generated from the frozen receipts
instead.  Every emitted file is written with ``open('x')`` so a rerun cannot
overwrite a submitted script.

This module writes shell scripts.  It reads no expression value and no outcome.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT_LITERAL = (
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/"
    "MASLD_Model_Benchmark"
)
RUNTIME_LITERAL = "${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738"
RUNTIME_SHA256 = "f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb"
SUMMARY_LITERAL = "${ROOT}/executions/model-data-088-21112137-summarization/summary"
PARTICIPANTS_LITERAL = "${ROOT}/executions/model-data-072-21080651-participants"
PARTICIPANTS_SHA256 = "6655f23b2092e163532113e74d75d3c2a93857ee870b4f0a94fdcb4fa7f9dbc0"
GATE_LITERAL = (
    "${ROOT}/config/evaluation/gse83452_baseline_nash_transfer_promotion_gate.json"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


HEADER = """#!/usr/bin/env bash
#SBATCH --job-name={job_name}
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --output={root}/executions/%x-%j.out
#SBATCH --error={root}/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS={cpus}
export OPENBLAS_NUM_THREADS={cpus}
export MKL_NUM_THREADS={cpus}

ROOT={root}
RUNTIME={runtime}
RUNTIME_ARTIFACTS_SHA256={runtime_sha}
PYTHON=${{RUNTIME}}/env/bin/python
"""

TRAP = """TARGET=${{ROOT}}/executions/{job_name}-${{SLURM_JOB_ID}}-{suffix}
STAGE=${{TARGET}}.staging
FAILED=${{TARGET}}.failed

preserve_failed_stage() {{
  if [[ -d "${{STAGE}}" && ! -e "${{FAILED}}" ]]; then
    mv "${{STAGE}}" "${{FAILED}}"
  fi
}}
trap preserve_failed_stage EXIT
if [[ -e "${{TARGET}}" || -e "${{STAGE}}" || -e "${{FAILED}}"{extra_guard} ]]; then
  printf 'Refusing to overwrite {label} artifacts.\\n' >&2
  exit 1
fi
install -d -m 0750 "${{STAGE}}"
cd "${{ROOT}}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${{ROOT}}/src:${{ROOT}}"
printf '%s  %s\\n' "${{RUNTIME_ARTIFACTS_SHA256}}" "${{RUNTIME}}/ARTIFACTS.json" \\
  | sha256sum --check --strict
bash -n slurm/{script_name}
"""


def register_script(*, arm: str, job_name: str, source_fit: str, source_fit_sha: str,
                    fit_contract: str, activation: str, script_name: str) -> str:
    contract = (
        "${ROOT}/config/evaluation/"
        f"gse83452_baseline_nash_transfer_prediction_{arm}_v1.json"
    )
    body = HEADER.format(
        job_name=job_name, cpus=1, mem="4G", time="00:20:00", root=ROOT_LITERAL,
        runtime=RUNTIME_LITERAL, runtime_sha=RUNTIME_SHA256,
    )
    body += f"ARM_ID={arm}\n"
    body += f"SOURCE_FIT={source_fit}\n"
    body += f"SOURCE_FIT_ARTIFACTS_SHA256={source_fit_sha}\n"
    body += f"FIT_CONTRACT=${{ROOT}}/{fit_contract}\n"
    body += f"ACTIVATION={activation}\n"
    body += f"SUMMARY={SUMMARY_LITERAL}\n"
    body += f"CONTRACT={contract}\n"
    body += TRAP.format(
        job_name=job_name, suffix=arm, label="the GSE83452 prediction contract",
        extra_guard=' || -e "${CONTRACT}"', script_name=script_name,
    )
    body += """printf '%s  %s\\n' "${SOURCE_FIT_ARTIFACTS_SHA256}" "${SOURCE_FIT}/ARTIFACTS.json" \\
  | sha256sum --check --strict
"${PYTHON}" - "${SOURCE_FIT}" "${ACTIVATION}" "${SUMMARY}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
print("input frozen-tree verification passed")
PY
"${PYTHON}" scripts/register_gse83452_baseline_nash_transfer_prediction_contract.py \\
  --benchmark-root "${ROOT}" \\
  --source-fit "${SOURCE_FIT}" \\
  --source-fit-contract "${FIT_CONTRACT}" \\
  --activation "${ACTIVATION}" \\
  --external-expression "${SUMMARY}" \\
  --baseline-participants 152 \\
  --output "${CONTRACT}" \\
  > "${STAGE}/register_stdout.json" 2> "${STAGE}/register_stderr.txt"
cp -p "${CONTRACT}" "${STAGE}/prediction_contract.json"
sha256sum "${CONTRACT}" > "${STAGE}/contract.sha256"
"""
    body += f"""sha256sum \\
  scripts/register_gse83452_baseline_nash_transfer_prediction_contract.py \\
  slurm/{script_name} \\
  "${{FIT_CONTRACT}}" \\
  "${{SOURCE_FIT}}/ARTIFACTS.json" \\
  "${{ACTIVATION}}/ARTIFACTS.json" \\
  "${{SUMMARY}}/ARTIFACTS.json" \\
  > "${{STAGE}}/source.sha256"
"""
    body += """"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
contract = json.loads((root / "prediction_contract.json").read_text(encoding="utf-8"))
if (
    contract.get("status") != "registered_prediction_only_pending"
    or contract.get("external_labels_read") is not False
    or contract.get("external_metrics_allowed") is not False
    or contract.get("timepoint") != "baseline"
    or "age_sex_logistic" not in contract.get("unfittable_model_ids", [])
    or not contract.get("fitted_model_ids")
    or contract["selected_source_model_id"] in contract["unfittable_model_ids"]
):
    raise SystemExit("registered prediction contract differs")
freeze_tree(root, {
    "artifact_class": "gse83452_baseline_nash_transfer_prediction_contract_registration",
    "arm_id": contract["arm_id"],
    "selected_source_model_id": contract["selected_source_model_id"],
    "fitted_models": len(contract["fitted_model_ids"]),
    "unfittable_models": len(contract["unfittable_model_ids"]),
    "external_labels_read": False,
    "status": "passed",
})
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
cat "${TARGET}/contract.sha256"
"""
    return body


def predict_script(*, arm: str, job_name: str, contract_sha: str, source_fit: str,
                   activation: str, script_name: str) -> str:
    contract = (
        "${ROOT}/config/evaluation/"
        f"gse83452_baseline_nash_transfer_prediction_{arm}_v1.json"
    )
    body = HEADER.format(
        job_name=job_name, cpus=4, mem="24G", time="02:00:00", root=ROOT_LITERAL,
        runtime=RUNTIME_LITERAL, runtime_sha=RUNTIME_SHA256,
    )
    body += f"ARM_ID={arm}\n"
    body += f"CONTRACT={contract}\n"
    body += f"CONTRACT_SHA256={contract_sha}\n"
    body += f"SOURCE_FIT={source_fit}\n"
    body += f"ACTIVATION={activation}\n"
    body += f"SUMMARY={SUMMARY_LITERAL}\n"
    body += TRAP.format(
        job_name=job_name, suffix=arm, label="GSE83452 prediction-only",
        extra_guard="", script_name=script_name,
    )
    body += """printf '%s  %s\\n' "${CONTRACT_SHA256}" "${CONTRACT}" | sha256sum --check --strict
# The predictor must never name the record table or the evaluator outcome.
if grep -n 'gse83452_records\\.tsv\\|evaluator_only\\|labels\\.tsv' \\
    scripts/predict_gse83452_baseline_nash_transfer.py "${CONTRACT}"; then
  printf 'Prediction-only source unexpectedly names the outcome table.\\n' >&2
  exit 1
fi
"${PYTHON}" -m unittest -v \\
  tests.unit.test_predict_gse83452_baseline_nash_transfer \\
  tests.unit.test_fit_gse135251_nash_transfer_source_models \\
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/predict_gse83452_baseline_nash_transfer.py \\
  --benchmark-root "${ROOT}" \\
  --contract "${CONTRACT}" \\
  --contract-sha256 "${CONTRACT_SHA256}" \\
  --output "${STAGE}/predictions" \\
  > "${STAGE}/prediction_stdout.json" 2> "${STAGE}/prediction_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
"""
    body += f"""sha256sum \\
  scripts/predict_gse83452_baseline_nash_transfer.py \\
  scripts/fit_gse135251_nash_transfer_source_models.py \\
  tests/unit/test_predict_gse83452_baseline_nash_transfer.py \\
  slurm/{script_name} \\
  "${{CONTRACT}}" \\
  "${{RUNTIME}}/ARTIFACTS.json" \\
  "${{SOURCE_FIT}}/ARTIFACTS.json" \\
  "${{ACTIVATION}}/ARTIFACTS.json" \\
  "${{SUMMARY}}/ARTIFACTS.json" \\
  > "${{STAGE}}/source.sha256"
"""
    body += """"${PYTHON}" - "${STAGE}" "${ARM_ID}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.contracts import PredictionBundle

root = Path(sys.argv[1]).resolve(strict=True)
arm_id = sys.argv[2]
predictions = root / "predictions"
receipt = json.loads((predictions / "prediction_receipt.json").read_text(encoding="utf-8"))
with (predictions / "prediction_bundle_index.tsv").open(encoding="utf-8", newline="") as handle:
    bundles = list(csv.DictReader(handle, delimiter="\\t"))
if (
    receipt.get("status") != "passed_frozen_prediction_only_no_label_access"
    or receipt.get("arm_id") != arm_id
    or receipt.get("participants") != 152
    or receipt.get("timepoint") != "baseline"
    or receipt.get("followup_arrays_read") is not False
    or receipt.get("shared_genes") != 21465
    or receipt.get("external_labels_read") is not False
    or receipt.get("external_participant_table_read") is not False
    or receipt.get("query_fit_or_calibration_performed") is not False
    or receipt.get("external_metrics_calculated") is not False
    or "age_sex_logistic" not in receipt.get("unfittable_model_ids", [])
    or len(bundles) != receipt["prediction_bundles"]
):
    raise SystemExit("GSE83452 prediction-only verification differs")
for row in bundles:
    bundle_root = predictions / row["model_id"]
    verify_frozen_tree(bundle_root)
    bundle = PredictionBundle.load_json(bundle_root / "prediction_bundle.json")
    bundle.validate_artifacts(bundle_root)
    if (
        bundle.n_predictions != 152
        or bundle.metadata.get("external_labels_read") is not False
        or bundle.metadata.get("query_fit_or_calibration_performed") is not False
        or bundle.metadata.get("prediction_frozen_before_evaluator_label_join") is not True
    ):
        raise SystemExit("GSE83452 PredictionBundle differs")
    with (bundle_root / "predictions.tsv").open(encoding="utf-8", newline="") as handle:
        table = list(csv.DictReader(handle, delimiter="\\t"))
    if {"nash_status", "age", "sex", "intervention"} & set(table[0]):
        raise SystemExit("prediction table carries an outcome or covariate column")
freeze_tree(
    predictions,
    {
        "artifact_class": "gse83452_baseline_nash_transfer_predictions",
        "arm_id": arm_id,
        "participants": 152,
        "prediction_bundles": len(bundles),
        "external_labels_read": False,
        "query_fit_or_calibration_performed": False,
        "metrics_calculated": False,
        "status": "passed_unscored",
    },
)
verify_frozen_tree(predictions)
freeze_tree(
    root,
    {
        "artifact_class": "gse83452_baseline_nash_transfer_prediction_campaign",
        "arm_id": arm_id,
        "external_labels_read": False,
        "query_fit_or_calibration_performed": False,
        "metrics_calculated": False,
        "status": "passed_unscored",
    },
)
verify_frozen_tree(root)
print(json.dumps({
    "models": receipt["models"],
    "distinct_seed_prediction_vectors": receipt["distinct_seed_prediction_vectors"],
    "representation_spread_diagnostics": receipt["representation_spread_diagnostics"],
}, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/predictions/ARTIFACTS.json"
for model in $(awk -F'\\t' 'NR>1{print $1}' "${TARGET}/predictions/prediction_bundle_index.tsv"); do
  sha256sum "${TARGET}/predictions/${model}/ARTIFACTS.json"
done
"""
    return body


def labels_script(*, arm: str, job_name: str, predictions: str, models: list[str],
                  shas: list[str], script_name: str) -> str:
    body = HEADER.format(
        job_name=job_name, cpus=1, mem="4G", time="00:20:00", root=ROOT_LITERAL,
        runtime=RUNTIME_LITERAL, runtime_sha=RUNTIME_SHA256,
    )
    body += f"ARM_ID={arm}\n"
    body += f"PARTICIPANTS={PARTICIPANTS_LITERAL}\n"
    body += f"PARTICIPANTS_ARTIFACTS_SHA256={PARTICIPANTS_SHA256}\n"
    body += f"PREDICTIONS={predictions}\n"
    body += TRAP.format(
        job_name=job_name, suffix=f"{arm}-labels",
        label="GSE83452 evaluator-only label", extra_guard="", script_name=script_name,
    )
    body += """printf '%s  %s\\n' "${PARTICIPANTS_ARTIFACTS_SHA256}" "${PARTICIPANTS}/ARTIFACTS.json" \\
  | sha256sum --check --strict
install -d -m 0750 "${STAGE}/evaluator_only"
"${PYTHON}" -m unittest -v \\
  tests.unit.test_build_gse83452_nash_evaluator_labels \\
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/build_gse83452_nash_evaluator_labels.py \\
  --records "${PARTICIPANTS}/gse83452_records.tsv" \\
"""
    for model, digest in zip(models, shas, strict=True):
        body += f'  --prediction "${{PREDICTIONS}}/{model}" \\\n'
        body += f"  --prediction-artifacts-sha256 {digest} \\\n"
    body += """  --output-labels "${STAGE}/evaluator_only/labels.tsv" \\
  --output-receipt "${STAGE}/evaluator_only/label_receipt.json" \\
  > "${STAGE}/label_stdout.json" 2> "${STAGE}/label_stderr.txt"
"""
    body += f"""sha256sum \\
  scripts/build_gse83452_nash_evaluator_labels.py \\
  tests/unit/test_build_gse83452_nash_evaluator_labels.py \\
  slurm/{script_name} \\
  "${{PARTICIPANTS}}/ARTIFACTS.json" \\
  > "${{STAGE}}/source.sha256"
"""
    body += """"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
evaluator = root / "evaluator_only"
receipt = json.loads((evaluator / "label_receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_evaluator_only_outcome_materialized_after_prediction_freeze"
    or receipt.get("baseline_records") != 152
    or receipt.get("endpoint_evaluable") != 148
    or receipt.get("class_counts") != {"no_nash": 44, "nash": 104}
    or receipt.get("undefined_retained_as_missing") != 4
    or receipt.get("undefined_mapped_to_no_nash") is not False
):
    raise SystemExit("GSE83452 evaluator-only label verification differs")
freeze_tree(evaluator, {
    "artifact_class": "gse83452_baseline_nash_evaluator_only_labels",
    "endpoint_evaluable": 148,
    "outcome_is_evaluator_only": True,
    "materialized_after_prediction_freeze": True,
    "status": "passed",
})
verify_frozen_tree(evaluator)
freeze_tree(root, {
    "artifact_class": "gse83452_baseline_nash_evaluator_label_campaign",
    "metrics_calculated": False,
    "status": "passed",
})
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/evaluator_only/ARTIFACTS.json"
"""
    return body


def evaluate_script(*, arm: str, job_name: str, gate_eligible: bool, candidate: str,
                    unfittable: list[str], evaluator_only: str, evaluator_sha: str,
                    gate_sha: str, predictions: str, models: list[str],
                    shas: list[str], script_name: str) -> str:
    body = HEADER.format(
        job_name=job_name, cpus=8, mem="24G", time="06:00:00", root=ROOT_LITERAL,
        runtime=RUNTIME_LITERAL, runtime_sha=RUNTIME_SHA256,
    )
    body += f"ARM_ID={arm}\n"
    body += f"GATE={GATE_LITERAL}\n"
    body += f"GATE_SHA256={gate_sha}\n"
    body += f"EVALUATOR_ONLY={evaluator_only}\n"
    body += f"EVALUATOR_ONLY_ARTIFACTS_SHA256={evaluator_sha}\n"
    body += f"PREDICTIONS={predictions}\n"
    body += TRAP.format(
        job_name=job_name, suffix=f"{arm}-eval", label="GSE83452 evaluation",
        extra_guard="", script_name=script_name,
    )
    body += """printf '%s  %s\\n' "${GATE_SHA256}" "${GATE}" | sha256sum --check --strict
"${PYTHON}" -m unittest -v \\
  tests.unit.test_evaluate_gse83452_baseline_nash_transfer \\
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/evaluate_gse83452_baseline_nash_transfer.py \\
  --evaluator-only "${EVALUATOR_ONLY}" \\
  --evaluator-only-artifacts-sha256 "${EVALUATOR_ONLY_ARTIFACTS_SHA256}" \\
  --gate "${GATE}" \\
  --gate-sha256 "${GATE_SHA256}" \\
"""
    body += f'  --candidate-model-id {candidate} \\\n'
    body += f'  --arm-id {arm} \\\n'
    body += f'  --gate-eligible-arm {"true" if gate_eligible else "false"} \\\n'
    for model_id in unfittable:
        body += f"  --unfittable-baseline {model_id} \\\n"
    for model, digest in zip(models, shas, strict=True):
        body += f'  --prediction "${{PREDICTIONS}}/{model}" \\\n'
        body += f"  --prediction-artifacts-sha256 {digest} \\\n"
    body += """  --output "${STAGE}/evaluation" \\
  > "${STAGE}/evaluation_stdout.json" 2> "${STAGE}/evaluation_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
"""
    body += f"""sha256sum \\
  scripts/evaluate_gse83452_baseline_nash_transfer.py \\
  tests/unit/test_evaluate_gse83452_baseline_nash_transfer.py \\
  slurm/{script_name} \\
  "${{GATE}}" \\
  "${{EVALUATOR_ONLY}}/ARTIFACTS.json" \\
  > "${{STAGE}}/source.sha256"
"""
    body += """"${PYTHON}" - "${STAGE}" "${ARM_ID}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
arm_id = sys.argv[2]
evaluation = root / "evaluation"
receipt = json.loads((evaluation / "receipt.json").read_text(encoding="utf-8"))
verdict = json.loads((evaluation / "promotion_gate_verdict.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "passed_external_development_scoring"
    or receipt.get("arm_id") != arm_id
    or receipt.get("participants") != 148
    or receipt.get("frozen_predictions") != 152
    or receipt.get("undefined_mapped_to_no_nash") is not False
    or receipt.get("prevalence_is_not_the_auprc_null") is not True
    or receipt.get("champion_claim_eligible") is not False
    or verdict.get("verdict") not in {"PASS", "FAIL"}
):
    raise SystemExit("GSE83452 evaluation verification differs")
freeze_tree(evaluation, {
    "artifact_class": "gse83452_baseline_nash_transfer_evaluation",
    "arm_id": arm_id,
    "participants": 148,
    "promotion_gate_verdict": verdict["verdict"],
    "external_development_only": True,
    "champion_eligible": False,
    "status": "passed",
})
verify_frozen_tree(evaluation)
freeze_tree(root, {
    "artifact_class": "gse83452_baseline_nash_transfer_evaluation_campaign",
    "arm_id": arm_id,
    "promotion_gate_verdict": verdict["verdict"],
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps(verdict, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/evaluation/ARTIFACTS.json"
cat "${TARGET}/evaluation/promotion_gate_verdict.json"
"""
    return body


BUILDERS = {
    "register": register_script,
    "predict": predict_script,
    "labels": labels_script,
    "evaluate": evaluate_script,
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", required=True, choices=sorted(BUILDERS))
    parser.add_argument("--spec", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    spec = json.loads(arguments.spec.read_text(encoding="utf-8"))
    spec["script_name"] = arguments.output.name
    body = BUILDERS[arguments.step](**spec)
    with arguments.output.open("x", encoding="utf-8") as handle:
        handle.write(body)
    arguments.output.chmod(0o750)
    print(json.dumps({"path": str(arguments.output), "sha256": sha256_file(arguments.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
