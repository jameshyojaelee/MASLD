#!/usr/bin/env bash
#SBATCH --job-name=model-check-287
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG="${ROOT}/config/evaluation/sequence_cnn_bundle_21076966_partial_development.json"
SOURCE="${ROOT}/config/evaluation/sequence_cnn_bundle_21076966_partial_sources.tsv"
AUDIT="${ROOT}/scripts/audit_sequence_cnn_bundle_partial_contract.py"
PREPARE="${ROOT}/scripts/prepare_sequence_evaluation_matrix.py"
EVALUATE="${ROOT}/scripts/evaluate_sequence_accessibility.py"
TARGET="${ROOT}/executions/model-check-287-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed partial evaluator: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite partial evaluator: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/evaluations"

for binding in \
  "b81ef972b582d451daca266474cfe7bff51f4acdb7e2f727577f8b2f3afe1060 ${CONFIG}" \
  "96785e97074679c3a1e8922d3c3b36d4aecd426580fe6e84eeccb9b2314a7852 ${SOURCE}" \
  "c1cbf2dc2358eb3e981390f1ed2d7190e87b30b1b68c6369b028a7eb2066ec23 ${AUDIT}" \
  "6fb5cd3a47c78113bcd2c44b4168d16d95c7e7e1a0329428203b58ccd8cb0301 ${ROOT}/tests/unit/test_audit_sequence_cnn_bundle_partial_contract.py" \
  "5b835e722619e232db0aa0731cfededaa4494a4e117b0e609cf816f082a82258 ${PREPARE}" \
  "e0b9b25857263a05d855c11a3cbd90dd457acb7876b9f695dde274a5a1252be4 ${ROOT}/tests/unit/test_prepare_sequence_evaluation_matrix.py" \
  "be1b0dac61123d4875aea210d304ae6c2b67605d74079d490fc21f4a07419824 ${EVALUATE}" \
  "2a0c36437c57ebe34adfe3ed5185293f4682bf65a5d823f28f21f31cfa856b19 ${ROOT}/tests/unit/test_evaluate_sequence_accessibility.py"; do
  read -r expected path <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_sequence_cnn_bundle_partial_contract \
  tests.unit.test_prepare_sequence_evaluation_matrix \
  > "${STAGE}/validation/plan-tests.stdout.txt" \
  2> "${STAGE}/validation/plan-tests.stderr.txt"
PYTHONPATH="src:." python "${AUDIT}" \
  --config "${CONFIG}" \
  --project-root "${ROOT}" \
  --output "${STAGE}/coverage_contract.json" \
  > "${STAGE}/validation/coverage-audit.stdout.json" \
  2> "${STAGE}/validation/coverage-audit.stderr.txt"
PYTHONPATH="src:." python "${PREPARE}" \
  --source-manifest "${SOURCE}" \
  --project-root "${ROOT}" \
  --output "${STAGE}/plan" \
  > "${STAGE}/validation/prepare.stdout.json" \
  2> "${STAGE}/validation/prepare.stderr.txt"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_evaluate_sequence_accessibility \
  > "${STAGE}/validation/evaluator-tests.stdout.txt" \
  2> "${STAGE}/validation/evaluator-tests.stderr.txt"

while IFS=$'\t' read -r dataset_id lineage_id split_id manifest manifest_sha \
  observed observed_sha baseline baseline_sha candidate_count candidate_ids; do
  [[ "${dataset_id}" == gse296875 ]] || { printf 'Dataset differs.\n' >&2; exit 1; }
  if [[ "${split_id}" =~ ^donor([0-4])_genomic([0-4])$ ]]; then
    donor_fold="${BASH_REMATCH[1]}"
    genomic_fold="${BASH_REMATCH[2]}"
  else
    printf 'Donor/genomic split differs.\n' >&2
    exit 1
  fi
  if [[ "${donor_fold}" != "${genomic_fold}" ]]; then
    printf 'Donor/genomic split differs.\n' >&2
    exit 1
  fi
  [[ "${candidate_count}" == 1 ]] || { printf 'Candidate count differs.\n' >&2; exit 1; }
  [[ "${candidate_ids}" == sequence_cnn_control__s20260824 ]] || {
    printf 'Candidate seed identity differs.\n' >&2
    exit 1
  }
  manifest_path="${STAGE}/plan/${manifest}"
  printf '%s  %s\n' "${manifest_sha}" "${manifest_path}" | sha256sum --check --strict
  output="${STAGE}/evaluations/${lineage_id}/${split_id}"
  install -d -m 0750 "$(dirname "${output}")"
  PYTHONPATH="src:." python "${EVALUATE}" \
    --prediction-manifest "${manifest_path}" \
    --observed-profiles "${observed}" \
    --observed-artifacts-sha256 "${observed_sha}" \
    --training-baseline "${baseline}" \
    --baseline-artifacts-sha256 "${baseline_sha}" \
    --dataset-id "${dataset_id}" \
    --split-id "${split_id}" \
    --lineage-id "${lineage_id}" \
    --evaluation-role valid \
    --output "${output}" \
    > "${STAGE}/validation/evaluate-${lineage_id}-${split_id}.stdout.json" \
    2> "${STAGE}/validation/evaluate-${lineage_id}-${split_id}.stderr.txt"
done < <(tail -n +2 "${STAGE}/plan/matrix.tsv")

PYTHONPATH="src:." python - "${STAGE}" <<'PY'
from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

root = Path(sys.argv[1]).resolve(strict=True)
coverage = json.loads((root / "coverage_contract.json").read_text(encoding="utf-8"))
with (root / "plan/matrix.tsv").open(encoding="utf-8", newline="") as handle:
    plan = list(csv.DictReader(handle, delimiter="\t"))
if (
    coverage.get("status") != "pass_complete_eligible_eight_unit_partial_coverage"
    or coverage.get("unit_count") != 8
    or len(plan) != 8
):
    raise SystemExit("coverage/evaluation plan differs")
expected_units = {
    (row["lineage_id"], row["split_id"], int(row["seed"]))
    for row in coverage["units"]
}
plan_units = {
    (row["lineage_id"], row["split_id"], 20260824)
    for row in plan
}
if plan_units != expected_units:
    raise SystemExit("actual receipt seed or coverage differs")
evaluations = []
candidate = "sequence_cnn_control__s20260824"
for row in plan:
    relative = (
        Path("evaluations") / row["lineage_id"] / row["split_id"] / "evaluation.json"
    )
    result = json.loads((root / relative).read_text(encoding="utf-8"))
    if (
        result.get("status") != "pass"
        or result.get("dataset_id") != "gse296875"
        or result.get("lineage_id") != row["lineage_id"]
        or result.get("split_id") != row["split_id"]
        or result.get("evaluation_role") != "valid"
        or result.get("models")
        != ["uniform_global", "training_pseudobulk_mean", candidate]
        or result.get("candidate_valid_ranking") != [candidate]
        or result.get("candidate_seeds") != {candidate: 20260824}
        or list(result.get("prediction_artifacts", {})) != [candidate]
        or result.get("test_outcomes_read") is not False
        or result.get("champion_claim_allowed") is not False
        or result.get("promotion_gate_evaluated") is not False
        or result.get("external_evaluation") is not False
        or result.get("seeds_used_as_biological_replicates") is not False
        or result.get("cells_used_as_independent_replicates") is not False
    ):
        raise SystemExit(f"singleton evaluation safety contract differs: {relative}")
    evaluations.append(
        {
            "lineage_id": row["lineage_id"],
            "split_id": row["split_id"],
            "seed": 20260824,
            "candidate_id": candidate,
            "strongest_matching_training_control": result[
                "strongest_training_only_baseline_selected_on_valid"
            ],
            "result": relative.as_posix(),
        }
    )
summary = {
    "schema_version": "masld-bench-sequence-cnn-eight-unit-partial-evaluation-v1",
    "status": "pass_partial_development_only",
    "dataset_id": "gse296875",
    "evaluation_role": "valid",
    "biological_unit": "donor",
    "genomic_unit": "held_contig_block",
    "units": 8,
    "independent_evaluations": 8,
    "actual_seeds": [20260824],
    "evaluations": evaluations,
    "one_candidate_per_evaluator_call": True,
    "matching_controls": ["uniform_global", "training_pseudobulk_mean"],
    "donor_safe_diagonal_splits_only": True,
    "coverage_filled_post_hoc": False,
    "partial_development_only": True,
    "cross_model_ranking_calculated": False,
    "cross_family_ranking_calculated": False,
    "cross_unit_metric_aggregation_calculated": False,
    "promotion_gate_evaluated": False,
    "test_outcomes_read": False,
    "sealed_data_read": False,
    "external_evaluation": False,
    "champion_claim_allowed": False,
}
(root / "partial_evaluation_summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY

python -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${CONFIG}" \
  "${SOURCE}" \
  "${AUDIT}" \
  "${PREPARE}" \
  "${EVALUATE}" \
  "${ROOT}/tests/unit/test_audit_sequence_cnn_bundle_partial_contract.py" \
  "${ROOT}/tests/unit/test_prepare_sequence_evaluation_matrix.py" \
  "${ROOT}/tests/unit/test_evaluate_sequence_accessibility.py" \
  "${ROOT}/slurm/evaluate_sequence_cnn_bundle_21076966_partial_cpu.sbatch" \
  > "${STAGE}/source.sha256"
micromamba deactivate
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads(
    (root / "partial_evaluation_summary.json").read_text(encoding="utf-8")
)
if (
    summary["status"] != "pass_partial_development_only"
    or summary["units"] != 8
    or summary["independent_evaluations"] != 8
    or not summary["one_candidate_per_evaluator_call"]
    or summary["coverage_filled_post_hoc"]
    or summary["cross_model_ranking_calculated"]
    or summary["cross_family_ranking_calculated"]
    or summary["cross_unit_metric_aggregation_calculated"]
    or summary["promotion_gate_evaluated"]
    or summary["test_outcomes_read"]
    or summary["sealed_data_read"]
    or summary["external_evaluation"]
    or summary["champion_claim_allowed"]
):
    raise SystemExit("partial development summary differs")
freeze_tree(
    root,
    {
        "artifact_class": "sequence_cnn_eight_unit_partial_development_evaluation",
        "dataset_id": "gse296875",
        "evaluation_role": "valid",
        "units": 8,
        "independent_evaluations": 8,
        "actual_seeds": [20260824],
        "biological_unit": "donor",
        "genomic_unit": "held_contig_block",
        "one_candidate_per_evaluator_call": True,
        "coverage_filled_post_hoc": False,
        "cross_model_ranking_calculated": False,
        "cross_family_ranking_calculated": False,
        "cross_unit_metric_aggregation_calculated": False,
        "promotion_gate_evaluated": False,
        "test_outcomes_read": False,
        "sealed_data_read": False,
        "external_evaluation": False,
        "champion_claim_allowed": False,
        "status": "passed_partial_development_only",
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Partial development evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
