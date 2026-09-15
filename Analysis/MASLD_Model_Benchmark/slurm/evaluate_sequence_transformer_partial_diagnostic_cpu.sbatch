#!/usr/bin/env bash
#SBATCH --job-name=data-eval-254
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
CONFIG="${ROOT}/config/evaluation/sequence_transformer_partial_diagnostic_21076977.toml"
SOURCE="${ROOT}/config/evaluation/sequence_transformer_partial_diagnostic_21076977.tsv"
TARGET="${ROOT}/executions/data-eval-254-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed partial diagnostic: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite partial diagnostic: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/evaluations"

for binding in \
  "de959983d484b4e7c07c7105392ba536ac125369ea4d157beac2ef48150bd22c ${CONFIG}" \
  "87eb37ef30a10ad2429fa3b2df122b485e03c09d8da98b79d0ea875c08737e13 ${SOURCE}" \
  "be1b0dac61123d4875aea210d304ae6c2b67605d74079d490fc21f4a07419824 ${ROOT}/scripts/evaluate_sequence_accessibility.py" \
  "5b835e722619e232db0aa0731cfededaa4494a4e117b0e609cf816f082a82258 ${ROOT}/scripts/prepare_sequence_evaluation_matrix.py" \
  "2a0c36437c57ebe34adfe3ed5185293f4682bf65a5d823f28f21f31cfa856b19 ${ROOT}/tests/unit/test_evaluate_sequence_accessibility.py" \
  "e0b9b25857263a05d855c11a3cbd90dd457acb7876b9f695dde274a5a1252be4 ${ROOT}/tests/unit/test_prepare_sequence_evaluation_matrix.py" \
  "5166c44bafb15c2503ae8942947e63f726604bf45bd337fc16a949dade8c98b8 ${ROOT}/executions/sequence_gpu_bundle_023_sequence_transformer_control-ready_ea978e9249-21076977/ARTIFACTS.json"; do
  read -r expected path <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHON311="$(command -v python)"
PYTHONPATH="src:." python - "${ROOT}" "${CONFIG}" "${SOURCE}" "${STAGE}" <<'PY'
from pathlib import Path
import csv
from hashlib import sha256
import json
import sys
import tomllib

from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
config_path = Path(sys.argv[2]).resolve(strict=True)
source_path = Path(sys.argv[3]).resolve(strict=True)
output = Path(sys.argv[4]).resolve(strict=True)
config = tomllib.loads(config_path.read_text(encoding="utf-8"))
expected_units = set(config["expected_units"])
if (
    config["schema_version"] != "masld-bench-sequence-partial-diagnostic-v1"
    or config["status"] != "prospective_development_only_matched_unit_diagnostic"
    or config["expected_unit_count"] != 8
    or len(expected_units) != 8
    or config["candidate_count_per_unit"] != 1
    or config["expected_candidate_id"]
    != "sequence_transformer_control__s20260824"
    or config["chrombpnet_included"]
    or config["bpnet_included"]
    or config["sequence_cnn_included"]
    or config["cross_model_ranking_authorized"]
    or config["cross_fold_ranking_authorized"]
    or config["family_ranking_authorized"]
    or config["selection_authorized"]
    or config["promotion_authorized"]
    or config["champion_claim_authorized"]
    or config["test_outcome_access_authorized"]
    or config["sealed_data_access_authorized"]
    or config["external_evaluation"]
    or config["partial_coverage_may_imply_family_performance"]
):
    raise SystemExit("partial diagnostic authorization differs")

def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

bundle = root / config["prediction_bundle"]["path"]
if digest(bundle / "ARTIFACTS.json") != config["prediction_bundle"]["artifacts_sha256"]:
    raise SystemExit("prediction bundle differs")
verify_frozen_tree(bundle)
observed_units = set()
artifact_hashes = {}
for receipt_path in sorted((bundle / "task_receipts").glob("*.json")):
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    artifact = Path(receipt["expected_artifact"]).resolve(strict=True)
    artifact_manifest = json.loads(
        (artifact / "ARTIFACTS.json").read_text(encoding="utf-8")
    )
    metadata = artifact_manifest["metadata"]
    unit = f'{metadata["lineage_id"]}:{metadata["split_id"]}'
    if (
        receipt["status"] != "passed"
        or receipt["model_id"] != config["model_id"]
        or receipt["seed"] != config["seed"]
        or metadata["model_id"] != config["model_id"]
        or metadata["seed"] != config["seed"]
        or metadata["lineage_id"] != receipt["lineage_id"]
        or metadata["status"] != "passed"
        or unit in observed_units
    ):
        raise SystemExit("prediction task receipt differs")
    if digest(artifact / "ARTIFACTS.json") != receipt["artifacts_sha256"]:
        raise SystemExit("prediction task artifact differs")
    verify_frozen_tree(artifact)
    observed_units.add(unit)
    artifact_hashes[unit] = receipt["artifacts_sha256"]
if observed_units != expected_units:
    raise SystemExit("prediction unit census differs")

input_hashes = {}
for authority in config["evaluator_inputs"]:
    tree = root / authority["path"]
    if digest(tree / "ARTIFACTS.json") != authority["artifacts_sha256"]:
        raise SystemExit("evaluator-input bundle differs")
    verify_frozen_tree(tree)
    input_hashes[authority["lineage_id"]] = authority["artifacts_sha256"]
if set(input_hashes) != {unit.split(":", 1)[0] for unit in expected_units}:
    raise SystemExit("evaluator-input lineage census differs")

with source_path.open(encoding="utf-8", newline="") as handle:
    source_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(source_rows) != 5 or sum(row["artifact_role"] == "prediction_bundle" for row in source_rows) != 1:
    raise SystemExit("source manifest census differs")
receipt = {
    "schema_version": "masld-bench-sequence-partial-diagnostic-preflight-v1",
    "status": "pass_prospective_contract_before_outcome_read",
    "prediction_bundle_artifacts_sha256": config["prediction_bundle"]["artifacts_sha256"],
    "prediction_artifacts": dict(sorted(artifact_hashes.items())),
    "evaluator_input_bundles": dict(sorted(input_hashes.items())),
    "expected_units": sorted(expected_units),
    "candidate_id": config["expected_candidate_id"],
    "chrombpnet_included": False,
    "cross_model_ranking_authorized": False,
    "outcomes_read": False,
}
(output / "validation/preflight.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY

PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_prepare_sequence_evaluation_matrix \
  > "${STAGE}/validation/plan-tests.stdout.txt" \
  2> "${STAGE}/validation/plan-tests.stderr.txt"
PYTHONPATH="src:." python scripts/prepare_sequence_evaluation_matrix.py \
  --source-manifest "${SOURCE}" \
  --project-root "${ROOT}" \
  --output "${STAGE}/plan" \
  > "${STAGE}/validation/prepare.stdout.json" \
  2> "${STAGE}/validation/prepare.stderr.txt"
PYTHONPATH="src:." python - "${CONFIG}" "${STAGE}/plan/matrix.tsv" <<'PY'
from pathlib import Path
import csv
import sys
import tomllib

config = tomllib.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
with Path(sys.argv[2]).open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
units = {f'{row["lineage_id"]}:{row["split_id"]}' for row in rows}
if (
    units != set(config["expected_units"])
    or len(rows) != config["expected_unit_count"]
    or any(row["candidate_count"] != "1" for row in rows)
    or any(row["candidate_ids"] != config["expected_candidate_id"] for row in rows)
):
    raise SystemExit("prepared matched-unit matrix differs")
PY

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_evaluate_sequence_accessibility \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

while IFS=$'\t' read -r dataset_id lineage_id split_id manifest manifest_sha \
  observed observed_sha baseline baseline_sha candidate_count candidate_ids; do
  [[ "${dataset_id}" == gse296875 ]] || { printf 'Dataset identity differs.\n' >&2; exit 1; }
  [[ "${candidate_count}" == 1 ]] || { printf 'Candidate count differs.\n' >&2; exit 1; }
  [[ "${candidate_ids}" == sequence_transformer_control__s20260824 ]] \
    || { printf 'Candidate identity differs.\n' >&2; exit 1; }
  manifest_path="${STAGE}/plan/${manifest}"
  printf '%s  %s\n' "${manifest_sha}" "${manifest_path}" | sha256sum --check --strict
  output="${STAGE}/evaluations/${lineage_id}/${split_id}"
  install -d -m 0750 "$(dirname "${output}")"
  PYTHONPATH="src:." python scripts/evaluate_sequence_accessibility.py \
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

PYTHONPATH="src:." "${PYTHON311}" - "${CONFIG}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
import tomllib

config = tomllib.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
root = Path(sys.argv[2]).resolve(strict=True)
evaluations = []
observed_units = set()
for path in sorted((root / "evaluations").glob("*/*/evaluation.json")):
    result = json.loads(path.read_text(encoding="utf-8"))
    unit = f'{result["lineage_id"]}:{result["split_id"]}'
    if (
        result["status"] != "pass"
        or result["dataset_id"] != config["dataset_id"]
        or result["evaluation_role"] != "valid"
        or result["biological_unit"] != "donor"
        or result["genomic_unit"] != "held_contig_block"
        or sorted(result["prediction_artifacts"]) != [config["expected_candidate_id"]]
        or result["candidate_valid_ranking"] != [config["expected_candidate_id"]]
        or set(result["models"])
        != {"uniform_global", "training_pseudobulk_mean", config["expected_candidate_id"]}
        or result["test_outcomes_read"]
        or result["external_evaluation"]
        or result["champion_claim_allowed"]
        or result["promotion_gate_evaluated"]
        or result["seeds_used_as_biological_replicates"]
        or set(result["summaries"]) != {"valid"}
        or unit in observed_units
    ):
        raise SystemExit(f"partial diagnostic result differs: {path}")
    observed_units.add(unit)
    evaluations.append(
        {
            "unit": unit,
            "result": str(path.relative_to(root)),
            "candidate_id": config["expected_candidate_id"],
            "matching_baselines": ["uniform_global", "training_pseudobulk_mean"],
        }
    )
if observed_units != set(config["expected_units"]):
    raise SystemExit("partial diagnostic output coverage differs")
summary = {
    "schema_version": "masld-bench-sequence-partial-diagnostic-summary-v1",
    "status": "pass_partial_development_diagnostic_only",
    "dataset_id": config["dataset_id"],
    "candidate_id": config["expected_candidate_id"],
    "evaluation_role": "valid",
    "evaluation_units": len(evaluations),
    "evaluations": evaluations,
    "matched_unit_comparison_only": True,
    "chrombpnet_included": False,
    "bpnet_included": False,
    "sequence_cnn_included": False,
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
PY

python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${CONFIG}" \
  "${SOURCE}" \
  scripts/evaluate_sequence_accessibility.py \
  scripts/prepare_sequence_evaluation_matrix.py \
  tests/unit/test_evaluate_sequence_accessibility.py \
  tests/unit/test_prepare_sequence_evaluation_matrix.py \
  slurm/evaluate_sequence_transformer_partial_diagnostic_cpu.sbatch \
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
summary = json.loads((root / "diagnostic_summary.json").read_text(encoding="utf-8"))
if (
    summary["status"] != "pass_partial_development_diagnostic_only"
    or summary["evaluation_units"] != 8
    or not summary["matched_unit_comparison_only"]
    or summary["chrombpnet_included"]
    or summary["cross_model_ranking_calculated"]
    or summary["cross_fold_ranking_calculated"]
    or summary["family_ranking_calculated"]
    or summary["selection_performed"]
    or summary["promotion_gate_evaluated"]
    or summary["test_outcomes_read"]
    or summary["sealed_data_read"]
    or summary["external_evaluation"]
    or summary["champion_claim_allowed"]
):
    raise SystemExit("terminal partial diagnostic contract differs")
freeze_tree(
    root,
    {
        "artifact_class": "sequence_transformer_partial_development_diagnostic",
        "dataset_id": "gse296875",
        "model_id": "sequence_transformer_control",
        "seed": 20260824,
        "evaluation_role": "valid",
        "evaluation_units": 8,
        "biological_unit": "donor",
        "genomic_unit": "held_contig_block",
        "matched_unit_comparison_only": True,
        "chrombpnet_included": False,
        "cross_model_ranking_calculated": False,
        "cross_fold_ranking_calculated": False,
        "family_ranking_calculated": False,
        "selection_performed": False,
        "promotion_gate_evaluated": False,
        "test_outcomes_read": False,
        "sealed_data_read": False,
        "external_evaluation": False,
        "champion_claim_allowed": False,
        "status": "passed_partial_diagnostic_only",
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Partial sequence diagnostic: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
