#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-meta-032
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-meta-032-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-meta-032-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
FOLD0=${ROOT}/executions/scbasset-donor0-genomic0-seed11-valid-evaluation-21065092
FOLD0_SHA=7ab89892e508c145f7b3bae8819af0f4695624294cfcd22f1200984cdcba8341
FOLD1TO4=${ROOT}/executions/scbasset-fold1to4-valid-crossfold-evaluation-21070600
FOLD1TO4_SHA=2fb574e16edc290b19ca3420fc5258a766c1bce5a7e98550449420e0ac5f992c
AGGREGATOR=${ROOT}/scripts/aggregate_scbasset_all5_valid.py
AGGREGATOR_SHA=4eb8dd88410fd7812fac75970354dd5669b83828693b793d01b57690987cd8df
TEST=${ROOT}/tests/unit/test_aggregate_scbasset_all5_valid.py
TEST_SHA=dad81f0c6e0c6e90060b111745fa3fa929fe4196ddc425b1cccbbd3266c7dd81
TARGET=${ROOT}/executions/scbasset-all5-valid-meta-aggregate-${SLURM_JOB_ID}
STAGE=${TARGET}.staging

preserve_failed_stage() {
  local failed=${TARGET}.failed
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scBasset all-fold aggregate: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scBasset all-fold aggregate: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
printf '%s  %s\n' \
  "${FOLD0_SHA}" "${FOLD0}/ARTIFACTS.json" \
  "${FOLD1TO4_SHA}" "${FOLD1TO4}/ARTIFACTS.json" \
  "${AGGREGATOR_SHA}" "${AGGREGATOR}" \
  "${TEST_SHA}" "${TEST}" \
  | sha256sum --check --strict
PYTHONPATH=${ROOT}/src python - "${FOLD0}" "${FOLD1TO4}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
PYTHONPATH=${ROOT}:${ROOT}/src python -m py_compile "${AGGREGATOR}" "${TEST}"
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_aggregate_scbasset_all5_valid \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH=${ROOT}:${ROOT}/src python "${AGGREGATOR}" \
  --fold0-root "${FOLD0}" \
  --fold0-artifacts-sha256 "${FOLD0_SHA}" \
  --fold1to4-root "${FOLD1TO4}" \
  --fold1to4-artifacts-sha256 "${FOLD1TO4_SHA}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/aggregate.stdout.json"
sha256sum "${AGGREGATOR}" "${TEST}" "$0" > "${STAGE}/source.sha256"

PYTHONPATH=${ROOT}/src python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "evaluation/evaluation.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass"
    or receipt["outer_folds"] != [0, 1, 2, 3, 4]
    or receipt["n_donors"] != 39
    or receipt["donors_by_outer_fold"] != {"0": 9, "1": 7, "2": 4, "3": 8, "4": 11}
    or not receipt["mixed_seed_meta_aggregate"]
    or receipt["fixed_seed_five_fold_cv"]
    or not receipt["base_evaluator_identical_across_sources"]
    or not receipt["frozen_metric_rows_only"]
    or receipt["outcome_files_opened_or_recomputed"]
    or receipt["test_atac_read"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("scBasset all-fold aggregate contract differs")
print(freeze_tree(root, {
    "artifact_class": "scbasset_all5_valid_meta_aggregate",
    "model_id": "scbasset",
    "dataset_id": "gse296875",
    "outer_folds": [0, 1, 2, 3, 4],
    "n_donors": 39,
    "mixed_seed_meta_aggregate": True,
    "fold0_evaluation_artifacts_sha256": "7ab89892e508c145f7b3bae8819af0f4695624294cfcd22f1200984cdcba8341",
    "fold1to4_evaluation_artifacts_sha256": "2fb574e16edc290b19ca3420fc5258a766c1bce5a7e98550449420e0ac5f992c",
    "outcomes_recomputed": False,
    "test_outcomes_read": False,
    "champion_claim_allowed": False,
    "status": "pass",
}))
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scBasset all-fold meta-aggregate: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
