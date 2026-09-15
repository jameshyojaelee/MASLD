#!/usr/bin/env bash
#SBATCH --job-name=preflight-model-work-203
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
FIXTURE="${ROOT}/executions/gse281364-enformer-fixture-21070094"
TARGET="${ROOT}/executions/model-work-203-preflight-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
module load singularity/4.2.2
IMAGE="${ROOT}/executions/environments/enformer-crested-tensorflow-2.17.0-21065799/tensorflow-tensorflow-2.17.0-gpu.sif"
singularity exec --cleanenv --containall --net --network none \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 --env PYTHONPATH="${ROOT}" \
  --bind "${ROOT}:${ROOT}:ro" --pwd "${ROOT}" "${IMAGE}" \
  python -m unittest -v \
  tests.unit.test_build_gse281364_enformer_fixture \
  tests.unit.test_enformer_predict_gse281364 \
  > "${STAGE}/validation/tests.stdout.txt" 2> "${STAGE}/validation/tests.stderr.txt"
singularity exec --cleanenv --containall --net --network none \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 --env PYTHONPATH="${ROOT}" \
  --bind "${ROOT}:${ROOT}:ro" --pwd "${ROOT}" "${IMAGE}" \
  python scripts/enformer_predict_gse281364.py --help \
  > "${STAGE}/validation/direct-script-help.stdout.txt" \
  2> "${STAGE}/validation/direct-script-help.stderr.txt"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  executions/environments/enformer-crested-tensorflow-2.17.0-21065799 \
  executions/enformer-crested-restricted-port-admission-21065579 \
  executions/enformer-crested-safe-keras-inspection-21065759 \
  executions/enformer-crested-track-crosswalk-21065831 \
  executions/enformer-crested-restricted-forward-21066179 \
  executions/sequence-fasta-21063808 "${FIXTURE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
python - slurm/run_model_work_203.sbatch <<'PY'
from pathlib import Path
import sys
text = Path(sys.argv[1]).read_text()
required = (
    "#SBATCH --job-name=model-work-203",
    "#SBATCH --partition=gpu",
    "#SBATCH --gres=gpu:l40s:1",
    "#SBATCH --qos=nslab",
    "--containall --net --network none",
    '"${ROOT}:${ROOT}:ro"',
    '"${STAGE}:${STAGE}:rw"',
    '"${CACHE}:${CACHE}:rw"',
    '--pwd "${ROOT}"',
    '--env PYTHONPATH="${ROOT}"',
    '"${FIXTURE_AUDIT}"',
    '"${PREFLIGHT}"',
    'find "${STAGE}/scratch" -mindepth 1 -maxdepth 1',
    "29dc3835c1d13a6c6bd93315b554075107cc409a9179b11ee0e00180f67cef56",
    "6c195deb31303d5c1664341308f84625c3d862038fafc8725c781ce33d1c0dac",
    "44d4befd0e7faca5bf4eea57e94c296147d51818b805d6a8a8a2189ef2f4b457",
    "rm --",
)
for value in required:
    if text.count(value) != 1:
        raise SystemExit(f"GPU wrapper routing differs: {value}")
for forbidden in (
    "--qos=innovation", "masld-enformer", "PLACEHOLDER_", "rm -r", "rm -f"
):
    if forbidden in text:
        raise SystemExit(f"GPU wrapper contains forbidden value: {forbidden}")
PY
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum scripts/build_gse281364_enformer_fixture.py \
  scripts/enformer_predict_gse281364.py \
  scripts/enformer_crested_probe_checkpoint.py \
  tests/unit/test_build_gse281364_enformer_fixture.py \
  tests/unit/test_enformer_predict_gse281364.py \
  slurm/build_gse281364_enformer_fixture.sbatch \
  slurm/preflight_model_work_203_cpu.sbatch \
  slurm/run_model_work_203.sbatch > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
freeze_tree(Path(sys.argv[1]).resolve(strict=True), {
    "artifact_class": "enformer_gse281364_gpu_bundle_preflight",
    "dataset_id": "gse281364",
    "model_id": "enformer_crested_restricted_port",
    "gpu_submitted": False,
    "outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
