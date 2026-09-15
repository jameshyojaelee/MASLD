#!/usr/bin/env bash
#SBATCH --job-name=seq-eval-matrix
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
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
SOURCE="${ROOT}/config/evaluation/sequence_accessibility_matrix_completed_20260824.tsv"
TARGET="${ROOT}/executions/sequence-accessibility-valid-matrix-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed sequence evaluation matrix: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite sequence evaluation matrix: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/evaluations"

for binding in \
  "be1b0dac61123d4875aea210d304ae6c2b67605d74079d490fc21f4a07419824 ${ROOT}/scripts/evaluate_sequence_accessibility.py" \
  "5b835e722619e232db0aa0731cfededaa4494a4e117b0e609cf816f082a82258 ${ROOT}/scripts/prepare_sequence_evaluation_matrix.py" \
  "2a0c36437c57ebe34adfe3ed5185293f4682bf65a5d823f28f21f31cfa856b19 ${ROOT}/tests/unit/test_evaluate_sequence_accessibility.py" \
  "e0b9b25857263a05d855c11a3cbd90dd457acb7876b9f695dde274a5a1252be4 ${ROOT}/tests/unit/test_prepare_sequence_evaluation_matrix.py" \
  "0af83eb441de68b1f102a775340042e6dcdceecf37fcaeeea2662c4103e29a8c ${SOURCE}"; do
  read -r expected path <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
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

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_evaluate_sequence_accessibility \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

while IFS=$'\t' read -r dataset_id lineage_id split_id manifest manifest_sha \
  observed observed_sha baseline baseline_sha candidate_count candidate_ids; do
  [[ "${dataset_id}" == gse296875 ]] || { printf 'Dataset identity differs.\n' >&2; exit 1; }
  [[ "${candidate_count}" =~ ^[1-9][0-9]*$ ]] || { printf 'Candidate count differs.\n' >&2; exit 1; }
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

PYTHONPATH="src:." python - "${STAGE}" <<'PY'
from collections import defaultdict
import csv
import json
from pathlib import Path
import sys

root = Path(sys.argv[1]).resolve(strict=True)
with (root / "plan/matrix.tsv").open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
coverage = defaultdict(list)
evaluations = []
for row in rows:
    path = root / "evaluations" / row["lineage_id"] / row["split_id"] / "evaluation.json"
    result = json.loads(path.read_text(encoding="utf-8"))
    if (
        result["status"] != "pass"
        or result["dataset_id"] != row["dataset_id"]
        or result["lineage_id"] != row["lineage_id"]
        or result["split_id"] != row["split_id"]
        or result["evaluation_role"] != "valid"
        or result["test_outcomes_read"]
        or result["champion_claim_allowed"]
        or result["promotion_gate_evaluated"]
        or result["external_evaluation"]
        or result["seeds_used_as_biological_replicates"]
    ):
        raise SystemExit(f"evaluation safety contract differs: {path}")
    candidates = sorted(result["prediction_artifacts"])
    if candidates != sorted(row["candidate_ids"].split(",")):
        raise SystemExit(f"evaluation candidate census differs: {path}")
    group = f'{row["lineage_id"]}:{row["split_id"]}'
    for candidate in candidates:
        coverage[candidate].append(group)
    evaluations.append(
        {
            "lineage_id": row["lineage_id"],
            "split_id": row["split_id"],
            "candidate_ids": candidates,
            "result": str(path.relative_to(root)),
        }
    )
summary = {
    "schema_version": "masld-bench-sequence-accessibility-valid-matrix-v1",
    "status": "pass",
    "dataset_id": "gse296875",
    "evaluation_role": "valid",
    "biological_unit": "donor",
    "genomic_unit": "held_contig_block",
    "evaluation_groups": len(evaluations),
    "evaluations": evaluations,
    "candidate_coverage": dict(sorted(coverage.items())),
    "donor_balanced_within_each_fold": True,
    "folds_or_seeds_used_as_biological_replicates": False,
    "cross_fold_candidate_ranking_calculated": False,
    "incomplete_coverage_used_for_comparison": False,
    "test_outcomes_read": False,
    "sealed_data_read": False,
    "champion_claim_allowed": False,
    "external_evaluation": False,
}
(root / "matrix_summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY

python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  scripts/evaluate_sequence_accessibility.py \
  scripts/prepare_sequence_evaluation_matrix.py \
  tests/unit/test_evaluate_sequence_accessibility.py \
  tests/unit/test_prepare_sequence_evaluation_matrix.py \
  slurm/evaluate_sequence_accessibility_matrix.sbatch \
  config/evaluation/sequence_accessibility_matrix_completed_20260824.tsv \
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
summary = json.loads((root / "matrix_summary.json").read_text(encoding="utf-8"))
freeze_tree(
    root,
    {
        "artifact_class": "sequence_accessibility_valid_development_matrix",
        "dataset_id": "gse296875",
        "evaluation_role": "valid",
        "evaluation_groups": summary["evaluation_groups"],
        "candidate_ids": sorted(summary["candidate_coverage"]),
        "biological_unit": "donor",
        "genomic_unit": "held_contig_block",
        "folds_or_seeds_used_as_biological_replicates": False,
        "cross_fold_candidate_ranking_calculated": False,
        "test_outcomes_read": False,
        "sealed_data_read": False,
        "champion_claim_allowed": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Sequence accessibility valid matrix: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
