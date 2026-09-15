#!/usr/bin/env bash
#SBATCH --job-name=data-eval-261
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
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG=${ROOT}/config/evaluation/bpnet_partial_diagnostic_21077605.toml
CONFIG_SHA256=eb9e2121c329a9bfa50cab3045173c8f51849f9fda9a32d3a32b0045b50a18aa
SOURCE=${ROOT}/config/evaluation/bpnet_partial_diagnostic_21077605.tsv
SOURCE_SHA256=6ff090a7e0b56a635712d351102fee88fd00d32730105190f95a7c2f934366b7
BUNDLE=${ROOT}/executions/sequence_gpu_bundle_004_bpnet-ready_71d22f6453-21077605
BUNDLE_SHA256=14af86ac22e67bc29f8a7b85ee3a049b63f875ba5423210cc8e539d29b9fe9e8
PREPARE=${ROOT}/scripts/prepare_sequence_evaluation_matrix.py
PREPARE_SHA256=5b835e722619e232db0aa0731cfededaa4494a4e117b0e609cf816f082a82258
EVALUATOR=${ROOT}/scripts/evaluate_sequence_accessibility.py
EVALUATOR_SHA256=be1b0dac61123d4875aea210d304ae6c2b67605d74079d490fc21f4a07419824
PREPARE_TEST=${ROOT}/tests/unit/test_prepare_sequence_evaluation_matrix.py
PREPARE_TEST_SHA256=e0b9b25857263a05d855c11a3cbd90dd457acb7876b9f695dde274a5a1252be4
EVALUATOR_TEST=${ROOT}/tests/unit/test_evaluate_sequence_accessibility.py
EVALUATOR_TEST_SHA256=2a0c36437c57ebe34adfe3ed5185293f4682bf65a5d823f28f21f31cfa856b19
FAILED_ATTEMPT_ERR=${ROOT}/executions/data-eval-260-21097467.err
FAILED_ATTEMPT_ERR_SHA256=15a8c2e932cbf06cad13cc3206eab3b446e4f199c37ccc63b5c2d6ce589cfa82
FAILED_ATTEMPT_WRAPPER_SHA256=5ecec9f1027732ae33153919568e86d85212783d0de74ac35dcec6702884d809
TARGET=${ROOT}/executions/data-eval-261-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  status=$?
  if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${status}"
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite BPNet partial diagnostic artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/evaluations"
cd "${ROOT}"
for pair in \
  "${CONFIG_SHA256} ${CONFIG}" \
  "${SOURCE_SHA256} ${SOURCE}" \
  "${BUNDLE_SHA256} ${BUNDLE}/ARTIFACTS.json" \
  "${PREPARE_SHA256} ${PREPARE}" \
  "${EVALUATOR_SHA256} ${EVALUATOR}" \
  "${PREPARE_TEST_SHA256} ${PREPARE_TEST}" \
  "${EVALUATOR_TEST_SHA256} ${EVALUATOR_TEST}" \
  "${FAILED_ATTEMPT_ERR_SHA256} ${FAILED_ATTEMPT_ERR}"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHON311=$(command -v python)
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON311}" - "${ROOT}" "${CONFIG}" "${SOURCE}" "${BUNDLE}" \
  "${STAGE}/validation/preflight.json" <<'PY'
from pathlib import Path
import csv
import json
import sys
import tomllib

from masld_bench.artifacts import verify_frozen_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
config = tomllib.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
source = Path(sys.argv[3]).resolve(strict=True)
bundle = Path(sys.argv[4]).resolve(strict=True)
output = Path(sys.argv[5])
expected = dict(config["expected_candidates"])
if (
    config["schema_version"] != "masld-bench-sequence-partial-diagnostic-v1"
    or config["status"] != "prospective_development_only_matched_unit_diagnostic"
    or config["dataset_id"] != "gse296875"
    or config["model_id"] != "bpnet"
    or config["candidate_count_per_unit"] != 1
    or config["expected_unit_count"] != 8
    or len(expected) != 8
    or config["chrombpnet_included"]
    or not config["bpnet_included"]
    or config["sequence_cnn_included"]
    or config["sequence_transformer_included"]
    or any(config[key] for key in (
        "cross_model_ranking_authorized", "cross_fold_ranking_authorized",
        "family_ranking_authorized", "selection_authorized",
        "promotion_authorized", "champion_claim_authorized",
        "test_outcome_access_authorized", "sealed_data_access_authorized",
        "external_evaluation", "partial_coverage_may_imply_family_performance",
    ))
):
    raise SystemExit("BPNet partial diagnostic authorization differs")
manifest = verify_frozen_tree(bundle)
if (
    sha256_file(bundle / "ARTIFACTS.json")
    != config["prediction_bundle"]["artifacts_sha256"]
    or manifest["metadata"].get("bundle_id")
    != "sequence_gpu_bundle_004_bpnet-ready_71d22f6453"
    or manifest["metadata"].get("logical_tasks") != 8
    or manifest["metadata"].get("terminal_status") != "passed"
):
    raise SystemExit("BPNet prediction bundle differs")
observed = {}
artifact_hashes = {}
for receipt_path in sorted((bundle / "task_receipts").glob("*.json")):
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    artifact = Path(receipt["expected_artifact"]).resolve(strict=True)
    artifact_manifest = verify_frozen_tree(artifact)
    metadata = artifact_manifest["metadata"]
    unit = f'{metadata["lineage_id"]}:{metadata["split_id"]}'
    candidate = f'bpnet__s{metadata["seed"]}'
    if (
        receipt.get("status") != "passed"
        or receipt.get("model_id") != "bpnet"
        or metadata.get("model_id") != "bpnet"
        or metadata.get("dataset_id") != "gse296875"
        or metadata.get("status") != "passed"
        or sha256_file(artifact / "ARTIFACTS.json") != receipt["artifacts_sha256"]
        or expected.get(unit) != candidate
        or unit in observed
    ):
        raise SystemExit(f"BPNet task receipt differs: {receipt_path}")
    observed[unit] = candidate
    artifact_hashes[unit] = receipt["artifacts_sha256"]
if observed != expected:
    raise SystemExit("BPNet prediction unit census differs")
with source.open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if len(rows) != 6 or sum(row["artifact_role"] == "prediction_bundle" for row in rows) != 1:
    raise SystemExit("BPNet source manifest census differs")
output.write_text(json.dumps({
    "schema_version": "masld-bench-bpnet-partial-preflight-v1",
    "status": "pass_prospective_contract_before_outcome_read",
    "expected_candidates": dict(sorted(expected.items())),
    "prediction_artifacts": dict(sorted(artifact_hashes.items())),
    "outcomes_read": False,
    "sealed_data_read": False,
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

"${PYTHON311}" -m unittest -v \
  tests.unit.test_prepare_sequence_evaluation_matrix \
  > "${STAGE}/validation/plan-tests.stdout.txt" \
  2> "${STAGE}/validation/plan-tests.stderr.txt"
"${PYTHON311}" "${PREPARE}" \
  --source-manifest "${SOURCE}" --project-root "${ROOT}" \
  --output "${STAGE}/plan" \
  > "${STAGE}/validation/prepare.stdout.json" \
  2> "${STAGE}/validation/prepare.stderr.txt"
"${PYTHON311}" - "${CONFIG}" "${STAGE}/plan/matrix.tsv" <<'PY'
from pathlib import Path
import csv
import sys
import tomllib

config = tomllib.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
with Path(sys.argv[2]).open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
observed = {
    f'{row["lineage_id"]}:{row["split_id"]}': row["candidate_ids"]
    for row in rows
}
if (
    len(rows) != config["expected_unit_count"]
    or observed != dict(config["expected_candidates"])
    or any(row["candidate_count"] != "1" for row in rows)
):
    raise SystemExit("prepared BPNet matched-unit matrix differs")
PY

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
PYTHONPATH="${ROOT}:${ROOT}/src" python -m unittest -v \
  tests.unit.test_evaluate_sequence_accessibility \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
while IFS=$'\t' read -r dataset_id lineage_id split_id manifest manifest_sha \
  observed observed_sha baseline baseline_sha candidate_count candidate_ids; do
  [[ ${dataset_id} == gse296875 && ${candidate_count} == 1 ]] \
    || { printf 'BPNet matrix identity differs.\n' >&2; exit 1; }
  manifest_path=${STAGE}/plan/${manifest}
  printf '%s  %s\n' "${manifest_sha}" "${manifest_path}" | sha256sum --check --strict
  output=${STAGE}/evaluations/${lineage_id}/${split_id}
  install -d -m 0750 "$(dirname "${output}")"
  PYTHONPATH="${ROOT}:${ROOT}/src" python "${EVALUATOR}" \
    --prediction-manifest "${manifest_path}" \
    --observed-profiles "${observed}" \
    --observed-artifacts-sha256 "${observed_sha}" \
    --training-baseline "${baseline}" \
    --baseline-artifacts-sha256 "${baseline_sha}" \
    --dataset-id "${dataset_id}" --split-id "${split_id}" \
    --lineage-id "${lineage_id}" --evaluation-role valid \
    --output "${output}" \
    > "${STAGE}/validation/evaluate-${lineage_id}-${split_id}.stdout.json" \
    2> "${STAGE}/validation/evaluate-${lineage_id}-${split_id}.stderr.txt"
done < <(tail -n +2 "${STAGE}/plan/matrix.tsv")
python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
micromamba deactivate
module purge

module load python/3.11.5-GCCcore-13.2.0
sha256sum \
  "${CONFIG}" "${SOURCE}" "${PREPARE}" "${EVALUATOR}" \
  "${PREPARE_TEST}" "${EVALUATOR_TEST}" \
  "${FAILED_ATTEMPT_ERR}" \
  slurm/evaluate_bpnet_partial_diagnostic_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}:${ROOT}/src" python - "${CONFIG}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
import tomllib

from masld_bench.artifacts import freeze_tree

config = tomllib.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
root = Path(sys.argv[2]).resolve(strict=True)
expected = dict(config["expected_candidates"])
evaluations = []
observed = {}
for path in sorted((root / "evaluations").glob("*/*/evaluation.json")):
    result = json.loads(path.read_text(encoding="utf-8"))
    unit = f'{result["lineage_id"]}:{result["split_id"]}'
    candidate = expected.get(unit)
    if (
        candidate is None
        or result["status"] != "pass"
        or result["dataset_id"] != "gse296875"
        or result["evaluation_role"] != "valid"
        or result["biological_unit"] != "donor"
        or result["genomic_unit"] != "held_contig_block"
        or sorted(result["prediction_artifacts"]) != [candidate]
        or result["candidate_valid_ranking"] != [candidate]
        or set(result["models"])
        != {"uniform_global", "training_pseudobulk_mean", candidate}
        or result["test_outcomes_read"]
        or result["external_evaluation"]
        or result["champion_claim_allowed"]
        or result["promotion_gate_evaluated"]
        or result["seeds_used_as_biological_replicates"]
        or set(result["summaries"]) != {"valid"}
        or unit in observed
    ):
        raise SystemExit(f"BPNet partial result differs: {path}")
    observed[unit] = candidate
    evaluations.append({
        "unit": unit,
        "candidate_id": candidate,
        "result": str(path.relative_to(root)),
        "matching_baselines": ["uniform_global", "training_pseudobulk_mean"],
    })
if observed != expected:
    raise SystemExit("BPNet partial evaluation coverage differs")
summary = {
    "schema_version": "masld-bench-bpnet-partial-diagnostic-summary-v1",
    "status": "pass_partial_development_diagnostic_only",
    "dataset_id": "gse296875",
    "model_id": "bpnet",
    "evaluation_role": "valid",
    "evaluation_units": len(evaluations),
    "evaluations": evaluations,
    "matched_unit_comparison_only": True,
    "cross_model_ranking_calculated": False,
    "cross_fold_ranking_calculated": False,
    "family_ranking_calculated": False,
    "selection_performed": False,
    "promotion_gate_evaluated": False,
    "test_outcomes_read": False,
    "sealed_data_read": False,
    "external_evaluation": False,
    "champion_claim_allowed": False,
}
(root / "diagnostic_summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
freeze_tree(root, {
    "artifact_class": "bpnet_partial_development_diagnostic",
    "dataset_id": "gse296875", "model_id": "bpnet",
    "supersedes_failed_job": "21097467",
    "superseded_failure": "incorrect_per_unit_frozen_tree_assertion_after_all_evaluations",
    "superseded_wrapper_sha256": "5ecec9f1027732ae33153919568e86d85212783d0de74ac35dcec6702884d809",
    "evaluation_role": "valid", "evaluation_units": 8,
    "biological_unit": "donor", "genomic_unit": "held_contig_block",
    "matched_unit_comparison_only": True,
    "cross_model_ranking_calculated": False,
    "family_ranking_calculated": False, "selection_performed": False,
    "promotion_gate_evaluated": False, "test_outcomes_read": False,
    "sealed_data_read": False, "external_evaluation": False,
    "champion_claim_allowed": False, "status": "passed_partial_diagnostic_only",
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
