#!/usr/bin/env bash
#SBATCH --job-name=data-eval-247
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-eval-247-%j.out
#SBATCH --error=executions/logs/data-eval-247-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
CAMPAIGN=${ROOT}/executions/model-train-246-21092247
CAMPAIGN_SHA256=9c7ed333bc55d41fc3a056230c55d85156a79f8240a702f938f6da6822087eb8
HISEQ_CAMPAIGN_SHA256=f141a8012214dc1c2f393226cd306475e2c917b22f85e90d3027e32c81a292de
NOVASEQ_CAMPAIGN_SHA256=7997b10de21ace826c5b9a658bd828c2d79efb718ba11d125a5612ec6b1a9eb5
FIXTURE=${ROOT}/executions/model-data-080-21082249-gse274114
FIXTURE_SHA256=b334ff4fe77903be82bb228fc479a28d92d74f4cec0b06ea98c3d6a6100c9300
LABELS=${FIXTURE}/evaluator_only/source_labels.tsv
LABELS_SHA256=99bc28c64bbfeb3f73c55e386a3f2516e1b4ac8ec9e9518402497bcd601381b9
FOLDS=${FIXTURE}/evaluator_only/participant_folds.tsv
FOLDS_SHA256=4ba7cf2e06db2bbf445c8b0aeec4f76b1b56b576bf805cdf62d8843719b8aea3
HISEQ_TASK=${ROOT}/config/evaluation/gse274114_hiseq_healthy_vs_hbv_task.toml
HISEQ_TASK_SHA256=2a05bf09ac7fccd4b52cac787f089a956260bf4e150f29a53f0e4843c247d113
NOVASEQ_TASK=${ROOT}/config/evaluation/gse274114_novaseq_mash_vs_mash_hbv_task.toml
NOVASEQ_TASK_SHA256=6150de64a148d5b1fb8f07bdac783f70931924211a4486302cbfadc5ce471ce5
EVALUATOR_SOURCE_SHA256=2c436cfd946326276043856a8b8a3b201b86ec27c5e51cccbf4f340bb708dcca
EVALUATOR_TEST_SHA256=c78a3d3e83ccec49aeaa4b760814b9c1a86f0fd62c890bccb23c3c5f88bab62a
TARGET=${ROOT}/executions/data-eval-247-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite independent GSE274114 evaluation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CAMPAIGN_SHA256}" "${CAMPAIGN}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${HISEQ_CAMPAIGN_SHA256}" \
  "${CAMPAIGN}/campaigns/gse274114_hiseq_healthy_vs_hbv/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${NOVASEQ_CAMPAIGN_SHA256}" \
  "${CAMPAIGN}/campaigns/gse274114_novaseq_mash_vs_mash_hbv/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FIXTURE_SHA256}" "${FIXTURE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${LABELS_SHA256}" "${LABELS}" | sha256sum --check --strict
printf '%s  %s\n' "${FOLDS_SHA256}" "${FOLDS}" | sha256sum --check --strict
printf '%s  %s\n' "${HISEQ_TASK_SHA256}" "${HISEQ_TASK}" \
  | sha256sum --check --strict
printf '%s  %s\n' "${NOVASEQ_TASK_SHA256}" "${NOVASEQ_TASK}" \
  | sha256sum --check --strict
printf '%s  %s\n' "${EVALUATOR_SOURCE_SHA256}" \
  "${ROOT}/scripts/evaluate_gse274114_bundled_baselines.py" \
  | sha256sum --check --strict
printf '%s  %s\n' "${EVALUATOR_TEST_SHA256}" \
  "${ROOT}/tests/unit/test_evaluate_gse274114_bundled_baselines.py" \
  | sha256sum --check --strict

"${PYTHON}" -m unittest -v \
  tests.unit.test_evaluate_gse274114_bundled_baselines \
  > "${STAGE}/unittest.stdout.txt" \
  2> "${STAGE}/unittest.stderr.txt"
"${PYTHON}" - <<'PY'
from pathlib import Path

for relative in (
    "scripts/evaluate_gse274114_bundled_baselines.py",
    "tests/unit/test_evaluate_gse274114_bundled_baselines.py",
):
    path = Path(relative)
    compile(path.read_text(encoding="utf-8"), path.as_posix(), "exec")
PY
bash -n slurm/evaluate_gse274114_bundled_baselines_cpu.sbatch

"${PYTHON}" scripts/evaluate_gse274114_bundled_baselines.py \
  --campaign-root "${CAMPAIGN}" \
  --campaign-artifacts-sha256 "${CAMPAIGN_SHA256}" \
  --labels "${LABELS}" \
  --labels-sha256 "${LABELS_SHA256}" \
  --folds "${FOLDS}" \
  --folds-sha256 "${FOLDS_SHA256}" \
  --hiseq-task-spec "${HISEQ_TASK}" \
  --novaseq-task-spec "${NOVASEQ_TASK}" \
  --bootstrap-replicates 10000 \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/evaluation.stdout.json" \
  2> "${STAGE}/evaluation.stderr.txt"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1])
manifest = verify_frozen_tree(root / "evaluation")
receipt = json.loads((root / "evaluation/receipt.json").read_text(encoding="utf-8"))
with (root / "evaluation/registered_metric_intervals.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    intervals = list(csv.DictReader(handle, delimiter="\t"))
tests = (root / "unittest.stderr.txt").read_text(encoding="utf-8")
if (
    "Ran 6 tests" not in tests
    or "OK" not in tests
    or "Warning" in tests
    or (root / "evaluation.stderr.txt").read_text(encoding="utf-8")
    or manifest["metadata"].get("artifact_class")
    != "gse274114_independent_development_evaluation"
    or receipt.get("status")
    != "development_only_platform_specific_evaluation_complete"
    or receipt.get("prediction_bundles_verified_before_labels_opened") is not True
    or receipt.get("prediction_bundle_count") != 120
    or len(receipt.get("tasks", ())) != 2
    or len(receipt.get("model_results", ())) != 8
    or receipt.get("registered_metrics")
    != [
        "participant_macro_f1",
        "participant_auprc",
        "participant_auroc",
        "participant_brier",
        "participant_log_loss",
    ]
    or receipt.get("biological_unit") != "participant"
    or receipt.get("bootstrap_replicates") != 10000
    or receipt.get("nominal_or_confirmatory_p_values_computed") is not False
    or receipt.get("global_four_class_scoring_run") is not False
    or receipt.get("global_ood_scoring_run") is not False
    or receipt.get("pooled_cross_platform_scoring_run") is not False
    or receipt.get("task_champion_selected") is not False
    or receipt.get("sealed_or_champion_claim_eligible") is not False
    or len(intervals) != 40
):
    raise SystemExit("independent evaluation contract differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse274114_independent_evaluation_execution",
        "status": receipt["status"],
        "prediction_bundle_count": 120,
        "task_count": 2,
        "model_count": 8,
        "registered_metric_rows": 40,
        "biological_unit": "participant",
        "pooled_cross_platform_scoring_run": False,
        "task_champion_selected": False,
        "sealed_or_champion_claim_eligible": False,
    },
)
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/evaluate_gse274114_bundled_baselines.py \
  tests/unit/test_evaluate_gse274114_bundled_baselines.py \
  config/evaluation/gse274114_hiseq_healthy_vs_hbv_task.toml \
  config/evaluation/gse274114_novaseq_mash_vs_mash_hbv_task.toml \
  slurm/evaluate_gse274114_bundled_baselines_cpu.sbatch \
  "${CAMPAIGN}/ARTIFACTS.json" \
  "${CAMPAIGN}/campaigns/gse274114_hiseq_healthy_vs_hbv/ARTIFACTS.json" \
  "${CAMPAIGN}/campaigns/gse274114_novaseq_mash_vs_mash_hbv/ARTIFACTS.json" \
  "${FIXTURE}/ARTIFACTS.json" \
  "${LABELS}" \
  "${FOLDS}" \
  "${STAGE}/evaluation/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/evaluation/ARTIFACTS.json"
