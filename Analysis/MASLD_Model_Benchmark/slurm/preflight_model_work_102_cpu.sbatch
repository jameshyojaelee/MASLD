#!/usr/bin/env bash
#SBATCH --job-name=model-work-102-check
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
WRAPPER="${ROOT}/slurm/run_model_work_102.sbatch"
TARGET="${ROOT}/executions/model-work-102-preflight-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite model-work-102 preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' \
  363de52e6ceda82d9908d3a24fed47a18a164613eb16126c24ffbef0e564a225 \
  "${WRAPPER}" | sha256sum --check --strict
bash -n "${WRAPPER}"
PREFLIGHT_ONLY=1 bash "${WRAPPER}" \
  > "${STAGE}/validation/wrapper.stdout.txt" \
  2> "${STAGE}/validation/wrapper.stderr.txt"
module purge
module load python/3.12.3-GCCcore-13.3.0
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
CONVERSION="${ROOT}/executions/environments/dna-lm-torch2.6-conversion-21064945"
RUNTIME_PROBE="${ROOT}/executions/mpralegnet-native-probe-21066276"
export PYTHONPATH="${RUNTIME_PROBE}/offline-guard:${CONVERSION}/site-packages:${TORCH}/site-packages:${BASE}/site-packages:${ROOT}"
python -m unittest -v \
  tests.unit.test_build_gse281364_mpralegnet_fixture \
  tests.unit.test_mpralegnet_predict_gse281364 \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1])
stdout = (root / "validation/wrapper.stdout.txt").read_text()
if "model-work-102 CPU preflight passed; GPU execution remains disabled" not in stdout:
    raise SystemExit("model-work-102 wrapper preflight differs")
(root / "validation/preflight_receipt.json").write_text(
    json.dumps(
        {
            "schema_version": "masld-bench-model-work-102-preflight-v1",
            "status": "passed",
            "gpu_used": False,
            "gpu_submission_performed": False,
            "queue_enabled": False,
            "wrapper_sha256": "363de52e6ceda82d9908d3a24fed47a18a164613eb16126c24ffbef0e564a225",
            "cpu_prediction_eligible": False,
            "required_prediction_device": "l40s",
            "outcomes_read": False,
            "reporter_counts_read": False,
            "standalone_champion_eligible": False,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n"
)
PY
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/mpralegnet_predict_gse281364.py" \
  "${ROOT}/tests/unit/test_mpralegnet_predict_gse281364.py" \
  "${ROOT}/slurm/run_model_work_102.sbatch" \
  "${ROOT}/slurm/preflight_model_work_102_cpu.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/preflight_receipt.json").read_text())
if (
    receipt["status"] != "passed"
    or receipt["gpu_used"]
    or receipt["gpu_submission_performed"]
    or receipt["queue_enabled"]
    or receipt["cpu_prediction_eligible"]
    or receipt["required_prediction_device"] != "l40s"
    or receipt["outcomes_read"]
    or receipt["reporter_counts_read"]
    or receipt["standalone_champion_eligible"]
):
    raise SystemExit("model-work-102 preflight receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "model_work_102_cpu_preflight",
        "model_id": "mpralegnet_hepg2_test1_val2",
        "dataset_id": "gse281364",
        "gpu_used": False,
        "gpu_submission_performed": False,
        "queue_enabled": False,
        "required_prediction_device": "l40s",
        "outcomes_read": False,
        "standalone_champion_eligible": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'model-work-102 CPU preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
