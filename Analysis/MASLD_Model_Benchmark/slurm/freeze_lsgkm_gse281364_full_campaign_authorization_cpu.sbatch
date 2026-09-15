#!/usr/bin/env bash
#SBATCH --job-name=model-check-613
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG="${ROOT}/config/lsgkm_gse281364_dinucleotide_full_campaign.json"
SCRIPT="${ROOT}/scripts/freeze_lsgkm_gse281364_full_campaign_authorization.py"
TEST="${ROOT}/tests/unit/test_freeze_lsgkm_gse281364_full_campaign_authorization.py"
TARGET="${ROOT}/executions/model-check-613-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite LS-GKM campaign authorization: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

printf '%s  %s\n' 79d7ca90cde98db5152f544ad25075a8dec0d7f421b9ce710c144880db8aa236 "${CONFIG}" | sha256sum --check --strict
printf '%s  %s\n' ef27670174ef8697c76064e28cea171dda75650ebc106653230b1dae4fa7466c "${SCRIPT}" | sha256sum --check --strict
printf '%s  %s\n' d26bf0f0bd5aee6f165f6722d0d56a4d0bcaa775a11e253f39ef577bc423e9e8 "${TEST}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}" python -m unittest -v \
  tests.unit.test_freeze_lsgkm_gse281364_full_campaign_authorization \
  > "${STAGE}/tests.stdout.txt" 2> "${STAGE}/tests.stderr.txt"
PYTHONPATH="${ROOT}" python "${SCRIPT}" \
  --project-root "${ROOT}" --config "${CONFIG}" \
  --output "${STAGE}/authorization" \
  > "${STAGE}/authorization.stdout.json" 2> "${STAGE}/authorization.stderr.txt"
python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
sha256sum "${SCRIPT}" "${TEST}" "${CONFIG}" \
  "${ROOT}/slurm/freeze_lsgkm_gse281364_full_campaign_authorization_cpu.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "authorization/full_campaign_authorization.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "authorized_remaining_24_outcome_blind_fits"
    or receipt["remaining_fit_count"] != 24
    or receipt["bundle_count"] != 5
    or len(receipt["authorized_run_fits"]) != 24
    or receipt["benchmark_metrics_calculated"]
    or any(receipt[field] for field in ("outcomes_read", "reporter_counts_read", "control_prediction_values_read", "sealed_assets_read"))
):
    raise SystemExit("LS-GKM full-campaign authorization receipt differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_exact_dinucleotide_full_campaign_authorization",
    "design_id": receipt["design_id"],
    "remaining_fit_count": receipt["remaining_fit_count"],
    "bundle_count": receipt["bundle_count"],
    "reused_fit": receipt["reused_fit"],
    "benchmark_metrics_calculated": False,
    "outcomes_read": False,
    "reporter_counts_read": False,
    "control_prediction_values_read": False,
    "sealed_assets_read": False,
    "status": receipt["status"],
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'LS-GKM full-campaign authorization: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
