#!/usr/bin/env bash
#SBATCH --job-name=model-train-246
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=executions/logs/model-train-246-%j.out
#SBATCH --error=executions/logs/model-train-246-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
FIREWALL_SOURCE_SHA256=d3f87d29f1e01dc8e61399c51e62f5a147bc29c3361f5c35896ff14bdaedb80a
CAMPAIGN_SOURCE_SHA256=b7bdb99fe9a4fbe9c2eff0e4716afaf091ac13addb58c1781b45c68d2c2a646c
FIREWALL_TEST_SHA256=83aaa65387465be9cdd03a82202889475d769e848ab1394139371c1fd4d41b3f
CAMPAIGN_TEST_SHA256=5cb2a5324a703abde20ff2267a3354ef18c08377dcb552044db3ee04f0f83fc8
CAMPAIGN_CONFIG_SHA256=0d4cc96674171fb681278499466341c36317c8effc38e46ab816f5267ffab41c
VALIDATION=${ROOT}/executions/data-check-245-21092223
VALIDATION_SHA256=3bb7069577847df6565b9e16d25c8f60854e4b6f0a63a569d8dc0afbd3203e13
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
TARGET=${ROOT}/executions/model-train-246-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite bundled GSE274114 v2 campaign.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/campaigns"
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
printf '%s  %s\n' "${FIREWALL_SOURCE_SHA256}" \
  "${ROOT}/scripts/gse274114_baseline_firewall.py" | sha256sum --check --strict
printf '%s  %s\n' "${CAMPAIGN_SOURCE_SHA256}" \
  "${ROOT}/scripts/run_gse274114_baseline_task_campaign.py" | sha256sum --check --strict
printf '%s  %s\n' "${FIREWALL_TEST_SHA256}" \
  "${ROOT}/tests/unit/test_gse274114_baseline_firewall.py" | sha256sum --check --strict
printf '%s  %s\n' "${CAMPAIGN_TEST_SHA256}" \
  "${ROOT}/tests/unit/test_run_gse274114_baseline_task_campaign.py" | sha256sum --check --strict
printf '%s  %s\n' "${CAMPAIGN_CONFIG_SHA256}" \
  "${ROOT}/config/campaigns/gse274114_within_platform_baselines_v2.toml" \
  | sha256sum --check --strict
sha256sum --check --strict "${VALIDATION}/source.sha256"

"${PYTHON}" -m unittest -v \
  tests.unit.test_run_gse274114_baseline_task_campaign \
  > "${STAGE}/campaign_unittest.stdout.txt" \
  2> "${STAGE}/campaign_unittest.stderr.txt"
"${PYTHON}" - <<'PY'
from pathlib import Path

for relative in (
    "scripts/gse274114_baseline_firewall.py",
    "scripts/run_gse274114_baseline_task_campaign.py",
    "tests/unit/test_gse274114_baseline_firewall.py",
    "tests/unit/test_run_gse274114_baseline_task_campaign.py",
):
    path = Path(relative)
    compile(path.read_text(encoding="utf-8"), path.as_posix(), "exec")
PY
bash -n slurm/run_gse274114_within_platform_baselines_v2_cpu.sbatch

for task_index in 0 1; do
  task_id=${TASK_IDS[${task_index}]}
  task_spec=${TASK_SPECS[${task_index}]}
  task_sha256=${TASK_SPEC_SHA256[${task_index}]}
  printf '%s  %s\n' "${task_sha256}" "${task_spec}" | sha256sum --check --strict
  "${PYTHON}" scripts/run_gse274114_baseline_task_campaign.py \
    --feature-root "${FEATURES}" \
    --feature-artifacts-sha256 "${FEATURES_SHA256}" \
    --labels "${LABELS}" \
    --folds "${FOLDS}" \
    --task-id "${task_id}" \
    --task-spec "${task_spec}" \
    --task-spec-sha256 "${task_sha256}" \
    --publication-root "${TARGET}/campaigns/${task_id}" \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --output "${STAGE}/campaigns/${task_id}" \
    > "${STAGE}/${task_id}.stdout.json" \
    2> "${STAGE}/${task_id}.stderr.txt"
done

"${PYTHON}" - "${STAGE}" "${TARGET}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1])
target = Path(sys.argv[2])
task_ids = (
    "gse274114_hiseq_healthy_vs_hbv",
    "gse274114_novaseq_mash_vs_mash_hbv",
)
tests = (root / "campaign_unittest.stderr.txt").read_text(encoding="utf-8")
if "Ran 4 tests" not in tests or "OK" not in tests:
    raise SystemExit("bundled campaign unit tests differ")
inventory = []
for task_id in task_ids:
    campaign = root / "campaigns" / task_id
    manifest = verify_frozen_tree(campaign)
    receipt = json.loads((campaign / "campaign_receipt.json").read_text(encoding="utf-8"))
    with (campaign / "prediction_inventory.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if (
        (root / f"{task_id}.stderr.txt").read_text(encoding="utf-8")
        or manifest["metadata"].get("artifact_class")
        != "gse274114_prediction_only_baseline_task_campaign"
        or receipt.get("status") != "passed_prediction_only_baseline_campaign"
        or receipt.get("task_id") != task_id
        or receipt.get("prediction_bundle_count") != 60
        or receipt.get("prediction_bundle_census")
        != {
            "training_class_prior": 5,
            "gene_rank_nearest_centroid": 5,
            "hvg_pca_elastic_net": 25,
            "hvg_pca_linear_svm": 25,
        }
        or receipt.get("held_fold_labels_in_views") is not False
        or receipt.get("participant_covariates_in_views") is not False
        or receipt.get("preprocessing_fit_inside_training_folds_only") is not True
        or receipt.get("prediction_outputs_contain_observed_labels") is not False
        or receipt.get("scoring_or_evaluation_run") is not False
        or receipt.get("sealed_or_champion_claim_eligible") is not False
        or len(rows) != 60
    ):
        raise SystemExit(f"prediction-only task contract differs: {task_id}")
    inventory.append(
        {
            "task_id": task_id,
            "campaign_root": (target / "campaigns" / task_id).as_posix(),
            "campaign_artifacts_sha256": sha256_file(campaign / "ARTIFACTS.json"),
            "prediction_bundle_count": 60,
        }
    )
write_json_exclusive(
    root / "bundled_campaign_receipt.json",
    {
        "schema_version": "masld-bench-gse274114-bundled-baseline-campaign-v2",
        "status": "passed_prediction_only_bundled_campaign",
        "task_campaigns": inventory,
        "task_count": 2,
        "prediction_bundle_count": 120,
        "scoring_or_evaluation_run": False,
        "prediction_outputs_contain_observed_labels": False,
        "global_four_class_task_run": False,
        "global_ood_task_run": False,
        "sealed_or_champion_claim_eligible": False,
        "prior_v1_attempts_reused": False,
    },
)
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/gse274114_baseline_firewall.py \
  scripts/run_gse274114_baseline_task_campaign.py \
  tests/unit/test_gse274114_baseline_firewall.py \
  tests/unit/test_run_gse274114_baseline_task_campaign.py \
  config/campaigns/gse274114_within_platform_baselines_v2.toml \
  "${TASK_SPECS[0]}" \
  "${TASK_SPECS[1]}" \
  slurm/run_gse274114_within_platform_baselines_v2_cpu.sbatch \
  "${VALIDATION}/ARTIFACTS.json" \
  "${FEATURE_STAGE}/ARTIFACTS.json" \
  "${FEATURES}/ARTIFACTS.json" \
  "${FIXTURE}/ARTIFACTS.json" \
  "${STAGE}/campaigns/${TASK_IDS[0]}/ARTIFACTS.json" \
  "${STAGE}/campaigns/${TASK_IDS[1]}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1])
receipt = json.loads((root / "bundled_campaign_receipt.json").read_text(encoding="utf-8"))
freeze_tree(
    root,
    {
        "artifact_class": "gse274114_bundled_baseline_execution_v2",
        "status": receipt["status"],
        "task_count": 2,
        "prediction_bundle_count": 120,
        "scoring_or_evaluation_run": False,
        "prediction_outputs_contain_observed_labels": False,
        "prohibited_contrast_run": False,
        "prior_v1_attempts_reused": False,
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum \
  "${TARGET}/ARTIFACTS.json" \
  "${TARGET}/campaigns/${TASK_IDS[0]}/ARTIFACTS.json" \
  "${TARGET}/campaigns/${TASK_IDS[1]}/ARTIFACTS.json"
