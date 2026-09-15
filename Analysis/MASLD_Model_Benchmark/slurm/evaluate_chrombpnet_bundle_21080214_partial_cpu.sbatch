#!/usr/bin/env bash
#SBATCH --job-name=model-check-285
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
CONFIG="${ROOT}/config/evaluation/chrombpnet_bundle_21080214_partial_development.json"
PREPARE="${ROOT}/scripts/prepare_chrombpnet_bundle_partial_evaluation.py"
EVALUATE="${ROOT}/scripts/evaluate_sequence_accessibility.py"
TARGET="${ROOT}/executions/model-check-285-${SLURM_JOB_ID}"
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
  "1ca00cd93df799048a6c03cb20e1136679638979fe02d4bc26a256bb2cf1453b ${CONFIG}" \
  "880b95dbd30bbe4c5697f2290efe0f4aa2e629815e72cf2e732f664b2337a42b ${PREPARE}" \
  "6770ddce7a6c068c83008b9df871c4a3fc22914a4e9df6de1bb95f90be8310f7 ${ROOT}/tests/unit/test_prepare_chrombpnet_bundle_partial_evaluation.py" \
  "be1b0dac61123d4875aea210d304ae6c2b67605d74079d490fc21f4a07419824 ${EVALUATE}" \
  "2a0c36437c57ebe34adfe3ed5185293f4682bf65a5d823f28f21f31cfa856b19 ${ROOT}/tests/unit/test_evaluate_sequence_accessibility.py"; do
  read -r expected path <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_prepare_chrombpnet_bundle_partial_evaluation \
  > "${STAGE}/validation/prepare-tests.stdout.txt" \
  2> "${STAGE}/validation/prepare-tests.stderr.txt"
PYTHONPATH="src:." python "${PREPARE}" \
  --config "${CONFIG}" \
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

while IFS=$'\t' read -r dataset_id lineage_id split_id outer_fold seed logical_work_id \
  view_id manifest manifest_sha prediction_root prediction_sha observed observed_sha \
  baseline baseline_sha; do
  [[ "${dataset_id}" == gse296875 ]] || { printf 'Dataset differs.\n' >&2; exit 1; }
  [[ "${split_id}" == "donor${outer_fold}_genomic${outer_fold}" ]] || {
    printf 'Donor/genomic split differs.\n' >&2
    exit 1
  }
  [[ "${logical_work_id}" == "seq__chrombpnet__${lineage_id}__d${outer_fold}g${outer_fold}__s${seed}" ]] || {
    printf 'Logical work identity differs.\n' >&2
    exit 1
  }
  [[ "${view_id}" == chrombpnet_full || "${view_id}" == chrombpnet_nobias ]] || {
    printf 'Native view differs.\n' >&2
    exit 1
  }
  manifest_path="${STAGE}/plan/${manifest}"
  printf '%s  %s\n' "${manifest_sha}" "${manifest_path}" | sha256sum --check --strict
  output="${STAGE}/evaluations/${lineage_id}/${split_id}/${view_id}"
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
    > "${STAGE}/validation/evaluate-${lineage_id}-${split_id}-${view_id}.stdout.json" \
    2> "${STAGE}/validation/evaluate-${lineage_id}-${split_id}-${view_id}.stderr.txt"
done < <(tail -n +2 "${STAGE}/plan/evaluation_plan.tsv")

PYTHONPATH="src:." python - "${STAGE}" <<'PY'
from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

root = Path(sys.argv[1]).resolve(strict=True)
with (root / "plan/evaluation_plan.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    plan = list(csv.DictReader(handle, delimiter="\t"))
if len(plan) != 16:
    raise SystemExit("partial evaluation count differs")
units = {
    (
        row["logical_work_id"],
        row["lineage_id"],
        row["split_id"],
        int(row["outer_fold"]),
        int(row["seed"]),
    )
    for row in plan
}
if len(units) != 8 or {row["view_id"] for row in plan} != {
    "chrombpnet_full",
    "chrombpnet_nobias",
}:
    raise SystemExit("partial unit/view census differs")
evaluations = []
for row in plan:
    candidate = f'{row["view_id"]}__s{row["seed"]}'
    relative = (
        Path("evaluations")
        / row["lineage_id"]
        / row["split_id"]
        / row["view_id"]
        / "evaluation.json"
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
        or result.get("candidate_seeds") != {candidate: int(row["seed"])}
        or list(result.get("prediction_artifacts", {})) != [candidate]
        or result.get("prediction_artifacts", {}).get(candidate)
        != row["prediction_artifacts_sha256"]
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
            "logical_work_id": row["logical_work_id"],
            "lineage_id": row["lineage_id"],
            "split_id": row["split_id"],
            "outer_fold": int(row["outer_fold"]),
            "seed": int(row["seed"]),
            "view_id": row["view_id"],
            "candidate_id": candidate,
            "strongest_matching_training_control": result[
                "strongest_training_only_baseline_selected_on_valid"
            ],
            "result": relative.as_posix(),
        }
    )
summary = {
    "schema_version": "masld-bench-chrombpnet-eight-unit-partial-evaluation-v1",
    "status": "pass_partial_development_only",
    "dataset_id": "gse296875",
    "evaluation_role": "valid",
    "biological_unit": "donor",
    "genomic_unit": "held_contig_block",
    "units": 8,
    "independent_evaluations": 16,
    "views": ["chrombpnet_full", "chrombpnet_nobias"],
    "evaluations": evaluations,
    "one_candidate_view_per_evaluator_call": True,
    "matching_controls": ["uniform_global", "training_pseudobulk_mean"],
    "donor_safe_diagonal_splits_only": True,
    "actual_artifact_seeds_preserved": True,
    "partial_development_only": True,
    "cross_view_ranking_calculated": False,
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
  "${PREPARE}" \
  "${EVALUATE}" \
  "${ROOT}/tests/unit/test_prepare_chrombpnet_bundle_partial_evaluation.py" \
  "${ROOT}/tests/unit/test_evaluate_sequence_accessibility.py" \
  "${ROOT}/slurm/evaluate_chrombpnet_bundle_21080214_partial_cpu.sbatch" \
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
    or summary["independent_evaluations"] != 16
    or not summary["one_candidate_view_per_evaluator_call"]
    or summary["cross_view_ranking_calculated"]
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
        "artifact_class": "chrombpnet_eight_unit_partial_development_evaluation",
        "dataset_id": "gse296875",
        "evaluation_role": "valid",
        "units": 8,
        "independent_evaluations": 16,
        "views": ["chrombpnet_full", "chrombpnet_nobias"],
        "biological_unit": "donor",
        "genomic_unit": "held_contig_block",
        "one_candidate_view_per_evaluator_call": True,
        "cross_view_ranking_calculated": False,
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
