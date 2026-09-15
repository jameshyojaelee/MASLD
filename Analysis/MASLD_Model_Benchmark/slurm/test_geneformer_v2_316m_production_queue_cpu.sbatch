#!/usr/bin/env bash
#SBATCH --job-name=model-check-075
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
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/decima/bin/python
QUEUE=${ROOT}/config/campaigns/gpu_bundle_queue/052-model-training-503.json
QUEUE_SHA256=508c93a7e49e776c33a9ad5695ba605b427f2552bc653dd691b4d63dd0c171e3
WRAPPER=${ROOT}/slurm/run_geneformer_v2_316m_frozen_screen_50000.sbatch
WRAPPER_SHA256=dd4f94365ef262322776298890dd5d59fb37f68aa96c0261a217fe39ce3998bb
ONE_BATCH=${ROOT}/executions/geneformer-v2-316m-one-batch-l40s-21080950
ONE_BATCH_SHA256=050263e5cf675fb3b421a751cf06b57f294a322d7128b1ff63a2f6f633764f68
SCREEN_PREFLIGHT=${ROOT}/executions/geneformer-v2-316m-frozen-screen-preflight-21078697
SCREEN_PREFLIGHT_SHA256=c6d9b1c6f8ad22b7f7fc15d4ba8d18c6c60fbd2681e8668cfd70871949dd4f94
TARGET=${ROOT}/executions/model-check-075-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failure() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${exit_code}"
}
trap preserve_failure EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Geneformer production-queue validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
printf '%s  %s\n' \
  "${QUEUE_SHA256}" "${QUEUE}" \
  "${WRAPPER_SHA256}" "${WRAPPER}" \
  "${ONE_BATCH_SHA256}" "${ONE_BATCH}/ARTIFACTS.json" \
  "${SCREEN_PREFLIGHT_SHA256}" "${SCREEN_PREFLIGHT}/ARTIFACTS.json" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"
bash -n "${WRAPPER}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_gpu_bundle_dispatcher \
  tests.unit.test_geneformer_v2_316m_frozen_screen_50000 \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${PYTHON}" - "${QUEUE}" "${ROOT}" > "${STAGE}/validation/queue.json" <<'PY'
from pathlib import Path
import json
import sys

from scripts.gpu_bundle_dispatcher import load_item

item = load_item(Path(sys.argv[1]), Path(sys.argv[2]).resolve(strict=True))
if (
    not item["_ready"]
    or item["bundle_id"] != "model-training-503"
    or item["priority"] != 52
    or item["logical_tasks"] != 1
    or item["exports"].get("ONE_BATCH_ARTIFACTS_SHA256")
    != "050263e5cf675fb3b421a751cf06b57f294a322d7128b1ff63a2f6f633764f68"
    or item["exports"].get("SCREEN_PREFLIGHT_ARTIFACTS_SHA256")
    != "c6d9b1c6f8ad22b7f7fc15d4ba8d18c6c60fbd2681e8668cfd70871949dd4f94"
):
    raise SystemExit("Geneformer production queue contract differs")
print(json.dumps({"bundle_id": item["bundle_id"], "ready": True}, sort_keys=True))
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/gpu_bundle_dispatcher.py \
  scripts/extract_geneformer_v2_316m_frozen_screen_50000.py \
  tests/unit/test_gpu_bundle_dispatcher.py \
  tests/unit/test_geneformer_v2_316m_frozen_screen_50000.py \
  slurm/run_geneformer_v2_316m_frozen_screen_50000.sbatch \
  slurm/test_geneformer_v2_316m_production_queue_cpu.sbatch \
  config/campaigns/gpu_bundle_queue/052-model-training-503.json \
  > "${STAGE}/validation/source.sha256"
"${PYTHON}" - "${STAGE}" "${QUEUE_SHA256}" "${WRAPPER_SHA256}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
tests = (stage / "validation/tests.stderr.txt").read_text(encoding="utf-8")
if "OK" not in tests:
    raise SystemExit("Geneformer production queue tests did not pass")
freeze_tree(
    stage,
    {
        "artifact_class": "geneformer_v2_316m_production_queue_validation",
        "bundle_id": "model-training-503",
        "queue_sha256": sys.argv[2],
        "wrapper_sha256": sys.argv[3],
        "gpu_job_submitted_by_validation": False,
        "metrics_calculated": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
