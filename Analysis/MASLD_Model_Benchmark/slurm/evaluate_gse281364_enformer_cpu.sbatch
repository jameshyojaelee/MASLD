#!/usr/bin/env bash
#SBATCH --job-name=masld-enformer-mpra-eval
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON="${PROJECT}/.agent/scanpy_env/bin/python"
: "${PREDICTIONS_ROOT:?PREDICTIONS_ROOT is required}"
: "${PREDICTIONS_SHA256:?PREDICTIONS_SHA256 is required}"
OUTCOMES="${ROOT}/executions/gse281364-replicate-outcomes-21066307"
OUTCOMES_SHA=cb42bd4037127a4e4d6a406d6f7b00c706bc93b631cdf93083ce994809e791db
FIXTURE="${ROOT}/executions/gse281364-enformer-fixture-21070094"
FIXTURE_SHA=0faf9d48511a76c282bfc412195deebda63503fb5685cd3fb2bfd77c72505231
PREFLIGHT="${ROOT}/executions/enformer-evaluator-preflight-21070703"
PREFLIGHT_SHA=83985940e89e48b630761a46b54adb16752f941327c58e6fadddc654e0d682bc
TARGET="${ROOT}/executions/gse281364-enformer-evaluation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed Enformer evaluation stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Enformer evaluation artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

printf '%s  %s\n' "${PREFLIGHT_SHA}" "${PREFLIGHT}/ARTIFACTS.json" \
  | sha256sum --check --strict
for pair in \
  "5d618f34c4a3ec41397bd33f7eca8cf5380705ce98c52a48bc75f3f2fa92bc68 ${ROOT}/scripts/evaluate_gse281364_enformer.py" \
  "436f8694237d120d217c00031ae4f8f42e99dc90f12590ecb9e4b5739930eda8 ${ROOT}/tests/unit/test_evaluate_gse281364_enformer.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${PREFLIGHT}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(root)
metadata = json.loads((root / "ARTIFACTS.json").read_text())["metadata"]
if (
    metadata["artifact_class"]
    != "gse281364_enformer_evaluator_outcome_free_preflight"
    or metadata["synthetic_tests"] != 7
    or not metadata["cli_help_passed"]
    or metadata["project_outcomes_read"]
    or metadata["prediction_artifact_read"]
    or metadata["sealed_outcomes_read"]
    or metadata["status"] != "passed"
):
    raise SystemExit("Enformer evaluator preflight contract differs")
PY
module purge

PYTHONPATH="${ROOT}" "${PYTHON}" -m unittest -v \
  tests.unit.test_evaluate_gse281364_enformer \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}" "${PYTHON}" \
  "${ROOT}/scripts/evaluate_gse281364_enformer.py" \
  --predictions "${PREDICTIONS_ROOT}" \
  --predictions-sha256 "${PREDICTIONS_SHA256}" \
  --outcomes "${OUTCOMES}" \
  --outcomes-sha256 "${OUTCOMES_SHA}" \
  --fixture "${FIXTURE}" \
  --fixture-sha256 "${FIXTURE_SHA}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluation.stdout.json" \
  2> "${STAGE}/validation/evaluation.stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/evaluate_gse281364_enformer.py" \
  "${ROOT}/tests/unit/test_evaluate_gse281364_enformer.py" \
  "${ROOT}/slurm/evaluate_gse281364_enformer_cpu.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${PREDICTIONS_SHA256}" <<'PY'
from pathlib import Path
import csv
import gzip
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
prediction_sha = sys.argv[2]
receipt = json.loads((root / "evaluation/receipt.json").read_text())
with (root / "evaluation/aggregate_metrics.tsv").open(newline="", encoding="utf-8") as handle:
    metrics = list(csv.DictReader(handle, delimiter="\t"))
with (root / "evaluation/outer_fold_metrics.tsv").open(newline="", encoding="utf-8") as handle:
    folds = list(csv.DictReader(handle, delimiter="\t"))
with gzip.open(
    root / "evaluation/matched_endpoints.tsv.gz", "rt", newline="", encoding="utf-8"
) as handle:
    matched = list(csv.DictReader(handle, delimiter="\t"))
if (
    receipt["status"] != "pass_descriptive_development_assay_transfer"
    or receipt["prediction_artifacts_sha256"] != prediction_sha
    or receipt["elements"] != 1033
    or receipt["outer_locus_sequence_groups"] != 1033
    or receipt["experimental_replicates_per_context"] != 4
    or receipt["experimental_replicates_used_as_independent_donors"]
    or receipt["donor_count"] != 0
    or receipt["models_ranked"]
    or receipt["confirmatory_inference"]
    or receipt["native_sonnet_parity_established"]
    or not receipt["restricted_comparator"]
    or receipt["external_evaluation"]
    or receipt["sealed_outcomes_read"]
    or receipt["champion_eligible"]
    or receipt["clinical_claim_supported"]
    or len(metrics) != 8
    or len(folds) != 40
    or len(matched) != 8264
):
    raise SystemExit("Enformer development evaluation contract differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse281364_enformer_development_evaluation",
        "dataset_id": "gse281364",
        "model_id": "enformer_crested_restricted_port",
        "evaluator_preflight_artifacts_sha256": "83985940e89e48b630761a46b54adb16752f941327c58e6fadddc654e0d682bc",
        "prediction_artifacts_sha256": prediction_sha,
        "outcome_artifacts_sha256": "cb42bd4037127a4e4d6a406d6f7b00c706bc93b631cdf93083ce994809e791db",
        "fixture_artifacts_sha256": "0faf9d48511a76c282bfc412195deebda63503fb5685cd3fb2bfd77c72505231",
        "outcome_role": "exposed_development_MPRA_only",
        "experimental_replicates_per_context": 4,
        "replicates_used_as_independent_donors": False,
        "donor_count": 0,
        "models_ranked": False,
        "restricted_comparator": True,
        "external_evaluation": False,
        "champion_eligible": False,
        "status": "passed",
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 Enformer development evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
