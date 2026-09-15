#!/usr/bin/env bash
#SBATCH --job-name=model-training-333
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

: "${PEAKVI_CROSS_COHORT_AUTHORITY:?dependency-dispatch authority is required}"
: "${PEAKVI_CROSS_COHORT_CONFIG:?campaign config is required}"
: "${PEAKVI_CROSS_COHORT_CONFIG_SHA256:?campaign config hash is required}"
if [[ "${PEAKVI_CROSS_COHORT_AUTHORITY}" != peakvi_cross_cohort_v1 ]]; then
  printf 'PeakVI cross-cohort dispatch authority differs\n' >&2
  exit 2
fi

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
EXPECTED_CONFIG="${ROOT}/config/campaigns/peakvi_gse296875_to_gse281367_fit_predict_20260825.json"
RUNTIME="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
RUNTIME_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
SITE_PACKAGES="${RUNTIME}/site-packages"
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_ARTIFACTS_SHA256=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
TORCH_SITE_PACKAGES="${TORCH}/site-packages"
TARGET="${ROOT}/executions/peakvi-cross-cohort-fit-predict-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ "${PEAKVI_CROSS_COHORT_CONFIG}" != "${EXPECTED_CONFIG}" ]]; then
  printf 'PeakVI cross-cohort config path differs\n' >&2
  exit 2
fi
printf '%s  %s\n' "${PEAKVI_CROSS_COHORT_CONFIG_SHA256}" "${EXPECTED_CONFIG}" \
  | sha256sum --check --strict
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite PeakVI cross-cohort bundle: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/commits"

for record in \
  "${RUNTIME_ARTIFACTS_SHA256} ${RUNTIME}/ARTIFACTS.json" \
  "${TORCH_ARTIFACTS_SHA256} ${TORCH}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" "${TORCH}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY

module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${TORCH_SITE_PACKAGES}:${SITE_PACKAGES}:${ROOT}:${ROOT}/src" python \
  "${ROOT}/scripts/peakvi_cross_cohort_fit_predict.py" \
  --root "${ROOT}" \
  --config "${EXPECTED_CONFIG}" \
  --output "${STAGE}/predictions" \
  > "${STAGE}/validation/fit_predict.stdout.json"
PYTHONPATH="${TORCH_SITE_PACKAGES}:${SITE_PACKAGES}" python -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${EXPECTED_CONFIG}" \
  "${ROOT}/scripts/peakvi_cross_cohort_fit_predict.py" \
  "${ROOT}/scripts/commit_gse281367_atac_transport_predictions.py" \
  "${ROOT}/tests/unit/test_peakvi_cross_cohort_fit_predict.py" \
  "${ROOT}/slurm/run_peakvi_cross_cohort_bundle.sbatch" \
  > "${STAGE}/validation/source.sha256"

PYTHONPATH="${TORCH_SITE_PACKAGES}:${SITE_PACKAGES}:${ROOT}/src:${ROOT}" python - \
  "${ROOT}" "${EXPECTED_CONFIG}" "${STAGE}" <<'PY'
from hashlib import sha256
from pathlib import Path
import json
import subprocess
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
config_path = Path(sys.argv[2]).resolve(strict=True)
stage = Path(sys.argv[3]).resolve(strict=True)
config = json.loads(config_path.read_text(encoding="utf-8"))
receipt = json.loads((stage / "predictions/receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "passed"
    or receipt.get("native_output_family") != "joint_representation"
    or receipt.get("projected_output_family") != "masked_accessibility_count"
    or len(receipt.get("rotations", [])) != 2
    or receipt.get("condition_labels_read") is not False
    or receipt.get("phenotype_values_read") is not False
    or receipt.get("development_outcomes_read") is not False
    or receipt.get("sealed_data_read") is not False
    or receipt.get("metrics_calculated") is not False
    or receipt.get("champion_claim_allowed") is not False
):
    raise SystemExit("PeakVI fit-predict receipt differs")
rotation_ids = ("valid_context_test_target", "test_context_valid_target")
commits = []
for rotation_id in rotation_ids:
    prediction = stage / "predictions" / rotation_id
    manifest = verify_frozen_tree(prediction)
    prediction_sha = sha256((prediction / "ARTIFACTS.json").read_bytes()).hexdigest()
    if (
        manifest.get("metadata", {}).get("artifact_class")
        != "peakvi_gse281367_blind_prediction"
        or manifest.get("metadata", {}).get("rotation_id") != rotation_id
        or manifest.get("metadata", {}).get("native_output_family")
        != "joint_representation"
        or manifest.get("metadata", {}).get("metrics_calculated") is not False
    ):
        raise SystemExit(f"PeakVI prediction artifact differs: {rotation_id}")
    commit_output = stage / "commits" / rotation_id
    command = [
        sys.executable,
        str(root / "scripts/commit_gse281367_atac_transport_predictions.py"),
        "--prediction-root",
        str(prediction),
        "--prediction-artifacts-sha256",
        prediction_sha,
        "--axis-root",
        str(root / config["exchange_axis"]["path"]),
        "--axis-artifacts-sha256",
        config["exchange_axis"]["artifacts_sha256"],
        "--output",
        str(commit_output),
    ]
    subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE)
    committed = verify_frozen_tree(commit_output)
    if (
        committed.get("metadata", {}).get("artifact_class")
        != "gse281367_atac_transport_prediction_commit"
        or committed.get("metadata", {}).get("development_outcomes_read") is not False
        or committed.get("metadata", {}).get("metrics_calculated") is not False
    ):
        raise SystemExit(f"PeakVI blind commit differs: {rotation_id}")
    commits.append(
        {
            "rotation_id": rotation_id,
            "prediction_artifacts_sha256": prediction_sha,
            "commit_artifacts_sha256": sha256(
                (commit_output / "ARTIFACTS.json").read_bytes()
            ).hexdigest(),
        }
    )
(stage / "commit_receipt.json").write_text(
    json.dumps(
        {
            "schema_version": "masld-bench-peakvi-cross-cohort-commit-receipt-v1",
            "status": "passed",
            "commits": commits,
            "condition_labels_read": False,
            "phenotype_values_read": False,
            "development_outcomes_read": False,
            "sealed_data_read": False,
            "metrics_calculated": False,
            "thresholds_or_calibration_changed": False,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
freeze_tree(
    stage,
    {
        "artifact_class": "peakvi_cross_cohort_fit_predict_and_blind_commit",
        "model_id": "peakvi_source_fit_fixed_window_decoder",
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse281367",
        "rotations": list(rotation_ids),
        "native_output_family": "joint_representation",
        "projected_output_family": "masked_accessibility_count",
        "condition_labels_read": False,
        "phenotype_values_read": False,
        "development_outcomes_read": False,
        "sealed_data_read": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(stage)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'PeakVI cross-cohort artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
