#!/usr/bin/env bash
#SBATCH --job-name=model-train-243
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-1%2
#SBATCH --output=executions/logs/model-train-243-%A_%a.out
#SBATCH --error=executions/logs/model-train-243-%A_%a.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
FIREWALL_SOURCE_SHA256=26445b56aef2df6b5822a7d59d1eae6779492256dfad544275e56e9dff9506f3
CAMPAIGN_SOURCE_SHA256=b7bdb99fe9a4fbe9c2eff0e4716afaf091ac13addb58c1781b45c68d2c2a646c
CAMPAIGN_TEST_SHA256=5cb2a5324a703abde20ff2267a3354ef18c08377dcb552044db3ee04f0f83fc8
CAMPAIGN_CONFIG_SHA256=5226aaa1507edf8b1e292157bbadc23df356506e2436bb78e350ba05dba7760d
VALIDATION=${ROOT}/executions/data-check-241-21092012
VALIDATION_SHA256=40d2d58f94344e86a31764a229d66090836f58046f5536a090637e1645eb3737
FEATURE_STAGE=${ROOT}/executions/data-stage-242-21092033-gse274114-unlabeled
FEATURE_STAGE_SHA256=da8f7eb0196fbe98ecafe6cb32672efafeba5118d4f32205a1db23758bf998a0
FEATURES=${FEATURE_STAGE}/features
FEATURES_SHA256=ec18f8528e0459a939c6483d073b2b5bdee3f3559ba52a136f386fac85b1a540
FIXTURE=${ROOT}/executions/model-data-080-21082249-gse274114
FIXTURE_SHA256=b334ff4fe77903be82bb228fc479a28d92d74f4cec0b06ea98c3d6a6100c9300
LABELS=${FIXTURE}/evaluator_only/source_labels.tsv
FOLDS=${FIXTURE}/evaluator_only/participant_folds.tsv
TASK_IDS=(gse274114_hiseq_healthy_vs_hbv gse274114_novaseq_mash_vs_mash_hbv)
TASK_SPECS=(
  "${ROOT}/config/evaluation/gse274114_hiseq_healthy_vs_hbv_task.toml"
  "${ROOT}/config/evaluation/gse274114_novaseq_mash_vs_mash_hbv_task.toml"
)
TASK_SPEC_SHA256=(
  2a05bf09ac7fccd4b52cac787f089a956260bf4e150f29a53f0e4843c247d113
  6150de64a148d5b1fb8f07bdac783f70931924211a4486302cbfadc5ce471ce5
)
TASK_ID=${TASK_IDS[${SLURM_ARRAY_TASK_ID}]}
TASK_SPEC=${TASK_SPECS[${SLURM_ARRAY_TASK_ID}]}
TASK_SHA256=${TASK_SPEC_SHA256[${SLURM_ARRAY_TASK_ID}]}
TARGET=${ROOT}/executions/model-train-243-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite GSE274114 baseline campaign task.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${VALIDATION_SHA256}" "${VALIDATION}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FEATURE_STAGE_SHA256}" "${FEATURE_STAGE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FEATURES_SHA256}" "${FEATURES}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FIXTURE_SHA256}" "${FIXTURE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${TASK_SHA256}" "${TASK_SPEC}" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FIREWALL_SOURCE_SHA256}" \
  "${ROOT}/scripts/gse274114_baseline_firewall.py" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CAMPAIGN_SOURCE_SHA256}" \
  "${ROOT}/scripts/run_gse274114_baseline_task_campaign.py" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CAMPAIGN_TEST_SHA256}" \
  "${ROOT}/tests/unit/test_run_gse274114_baseline_task_campaign.py" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CAMPAIGN_CONFIG_SHA256}" \
  "${ROOT}/config/campaigns/gse274114_within_platform_baselines_v1.toml" \
  | sha256sum --check --strict
sha256sum --check --strict "${VALIDATION}/source.sha256"

"${PYTHON}" -m unittest -v \
  tests.unit.test_run_gse274114_baseline_task_campaign \
  > "${STAGE}/campaign_unittest.stdout.txt" \
  2> "${STAGE}/campaign_unittest.stderr.txt"
"${PYTHON}" - <<'PY'
from pathlib import Path

for relative in (
    "scripts/run_gse274114_baseline_task_campaign.py",
    "tests/unit/test_run_gse274114_baseline_task_campaign.py",
):
    path = Path(relative)
    compile(path.read_text(encoding="utf-8"), path.as_posix(), "exec")
PY
bash -n slurm/run_gse274114_within_platform_baselines_cpu.sbatch

"${PYTHON}" scripts/run_gse274114_baseline_task_campaign.py \
  --feature-root "${FEATURES}" \
  --feature-artifacts-sha256 "${FEATURES_SHA256}" \
  --labels "${LABELS}" \
  --folds "${FOLDS}" \
  --task-id "${TASK_ID}" \
  --task-spec "${TASK_SPEC}" \
  --task-spec-sha256 "${TASK_SHA256}" \
  --publication-root "${TARGET}/campaign" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --output "${STAGE}/campaign" \
  > "${STAGE}/campaign.stdout.json" \
  2> "${STAGE}/campaign.stderr.txt"
"${PYTHON}" - "${STAGE}" "${TASK_ID}" "${TARGET}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1])
task_id = sys.argv[2]
target = Path(sys.argv[3])
campaign = root / "campaign"
manifest = verify_frozen_tree(campaign)
receipt = json.loads((campaign / "campaign_receipt.json").read_text(encoding="utf-8"))
with (campaign / "prediction_inventory.tsv").open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
tests = (root / "campaign_unittest.stderr.txt").read_text(encoding="utf-8")
campaign_stderr = (root / "campaign.stderr.txt").read_text(encoding="utf-8")
expected_census = {
    "training_class_prior": 5,
    "gene_rank_nearest_centroid": 5,
    "hvg_pca_elastic_net": 25,
    "hvg_pca_linear_svm": 25,
}
if (
    "Ran 4 tests" not in tests
    or "OK" not in tests
    or campaign_stderr
    or manifest["metadata"].get("artifact_class")
    != "gse274114_prediction_only_baseline_task_campaign"
    or receipt.get("status") != "passed_prediction_only_baseline_campaign"
    or receipt.get("task_id") != task_id
    or receipt.get("prediction_bundle_count") != 60
    or receipt.get("prediction_bundle_census") != expected_census
    or receipt.get("fixed_stochastic_seeds") != [1103, 2207, 3301, 4409, 5501]
    or receipt.get("held_fold_labels_in_views") is not False
    or receipt.get("participant_covariates_in_views") is not False
    or receipt.get("cohort_wide_normalization_run") is not False
    or receipt.get("preprocessing_fit_inside_training_folds_only") is not True
    or receipt.get("prediction_outputs_contain_observed_labels") is not False
    or receipt.get("scoring_or_evaluation_run") is not False
    or receipt.get("global_four_class_task_run") is not False
    or receipt.get("global_ood_task_run") is not False
    or receipt.get("sealed_or_champion_claim_eligible") is not False
    or len(rows) != 60
):
    raise SystemExit("prediction-only baseline campaign contract differs")
for row in rows:
    published_root = Path(row["bundle_root"])
    try:
        relative = published_root.relative_to(target / "campaign")
    except ValueError as error:
        raise SystemExit("published prediction path escapes target") from error
    bundle_root = campaign / relative
    if row["artifacts_sha256"] != __import__(
        "masld_bench.hashing", fromlist=["sha256_file"]
    ).sha256_file(bundle_root / "ARTIFACTS.json"):
        raise SystemExit("prediction bundle identity differs")
    bundle_manifest = verify_frozen_tree(bundle_root)
    if (
        bundle_manifest["metadata"].get("query_labels_opened") is not False
        or bundle_manifest["metadata"].get("scoring_run") is not False
    ):
        raise SystemExit("prediction bundle firewall differs")
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/gse274114_baseline_firewall.py \
  scripts/run_gse274114_baseline_task_campaign.py \
  tests/unit/test_gse274114_baseline_firewall.py \
  tests/unit/test_run_gse274114_baseline_task_campaign.py \
  config/campaigns/gse274114_within_platform_baselines_v1.toml \
  "${TASK_SPEC}" \
  slurm/run_gse274114_within_platform_baselines_cpu.sbatch \
  "${VALIDATION}/ARTIFACTS.json" \
  "${FEATURE_STAGE}/ARTIFACTS.json" \
  "${FEATURES}/ARTIFACTS.json" \
  "${FIXTURE}/ARTIFACTS.json" \
  "${STAGE}/campaign/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${TASK_ID}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1])
task_id = sys.argv[2]
verify_frozen_tree(root / "campaign")
receipt = json.loads((root / "campaign/campaign_receipt.json").read_text(encoding="utf-8"))
freeze_tree(
    root,
    {
        "artifact_class": "gse274114_baseline_task_execution",
        "status": receipt["status"],
        "task_id": task_id,
        "prediction_bundle_count": 60,
        "outer_folds": 5,
        "scoring_or_evaluation_run": False,
        "held_fold_labels_in_views": False,
        "participant_covariates_in_views": False,
        "prohibited_contrast_run": False,
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/campaign/ARTIFACTS.json"
