#!/usr/bin/env bash
#SBATCH --job-name=masld-multivi-peakvi-probe
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
RUNTIME_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
SITE_PACKAGES="${RUNTIME}/site-packages"
TORCH_OVERLAY="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_OVERLAY_ARTIFACTS_SHA256=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
TORCH_SITE_PACKAGES="${TORCH_OVERLAY}/site-packages"
TARGET="${ROOT}/executions/multivi-peakvi-no-outcome-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then
      printf 'Cannot preserve failed MultiVI/PeakVI stage: %s\n' "${failed}" >&2
      return
    fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed MultiVI/PeakVI stage: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite MultiVI/PeakVI probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${TORCH_OVERLAY_ARTIFACTS_SHA256}" \
  "${TORCH_OVERLAY}/ARTIFACTS.json" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" "${TORCH_OVERLAY}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
verify_frozen_tree(Path(sys.argv[2]))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${TORCH_SITE_PACKAGES}:${SITE_PACKAGES}" python \
  "${ROOT}/scripts/multivi_peakvi_no_outcome_probe.py" \
  --output "${STAGE}/probe" \
  > "${STAGE}/validation/probe.stdout.json"
PYTHONPATH="${TORCH_SITE_PACKAGES}:${SITE_PACKAGES}" python -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/multivi_peakvi_no_outcome_probe.py" \
  "${ROOT}/slurm/multivi_peakvi_run_no_outcome_probe.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${STAGE}" "${RUNTIME_ARTIFACTS_SHA256}" \
  "${TORCH_OVERLAY_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "probe/receipt.json").read_text())
multivi = receipt["multivi"]
peakvi = receipt["peakvi"]
if (
    receipt["status"] != "pass"
    or receipt["scvi_tools"] != "1.5.0.post1"
    or receipt["torch"] != "2.6.0+cu124"
    or receipt["torch_cuda_build"] != "12.4"
    or "L40S" not in receipt["gpu"]
    or receipt["project_data_read"]
    or receipt["outcomes_read"]
    or receipt["labels_read"]
    or receipt["same_nucleus_project_fit_executed"]
    or receipt["gse244832_false_pairing"]
    or receipt["gse281367_false_pairing"]
    or receipt["champion_claim_allowed"]
    or not multivi["explicit_state_bridge_matches_upstream_row_sum_masks"]
    or not multivi["observed_all_zero_rows_rejected"]
    or multivi["query_atac_state"] != "structurally_missing"
    or multivi["query_atac_nonzero_values"] != 0
    or not multivi["query_composition_invariant"]
    or not multivi["repeat_forward_bit_identical"]
    or not multivi["backward_executed"]
    or not multivi["state_dict_resume_executed"]
    or not multivi["rna_conditioned_atac_runtime_candidate"]
    or not peakvi["observed_atac_required_at_query"]
    or peakvi["rna_conditioned_atac_eligible"]
    or peakvi["sealed_rna_only_eligible"]
    or not peakvi["repeat_forward_bit_identical"]
    or not peakvi["backward_executed"]
    or not peakvi["state_dict_resume_executed"]
):
    raise SystemExit("MultiVI/PeakVI no-outcome runtime contract differs")
freeze_tree(
    root,
    {
        "artifact_class": "multivi_peakvi_no_outcome_runtime_probe",
        "model_roster": ["multivi", "peakvi"],
        "runtime_artifacts_sha256": sys.argv[2],
        "torch_overlay_artifacts_sha256": sys.argv[3],
        "project_data_read": False,
        "outcomes_read": False,
        "multivi_explicit_state_bridge_passed": True,
        "multivi_project_fit_pending": True,
        "peakvi_observed_atac_development_only": True,
        "champion_claim_allowed": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'MultiVI/PeakVI no-outcome probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
