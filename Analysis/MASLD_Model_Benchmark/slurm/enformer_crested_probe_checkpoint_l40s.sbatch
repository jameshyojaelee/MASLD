#!/usr/bin/env bash
#SBATCH --job-name=masld-enformer-port-forward
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/enformer-crested-tensorflow-2.17.0-21065799"
RUNTIME_SHA256=3534aee9e54d8fe6dc09c01c63aaead911590b97a08a7a686d56024520885ec6
ADMISSION="${ROOT}/executions/enformer-crested-restricted-port-admission-21065579"
ADMISSION_SHA256=e69877c43f47ce7d40962129dae7819a66123ac3f5bef7e0283ac654f77006a5
INSPECTION="${ROOT}/executions/enformer-crested-safe-keras-inspection-21065759"
INSPECTION_SHA256=155537e853b730d440244581fda83ddff9ba8c13c9dabaa1f57bf841bb69196f
CROSSWALK="${ROOT}/executions/enformer-crested-track-crosswalk-21065831"
CROSSWALK_SHA256=a658a159b03fb05c9a468b4fdbc9c4d745916ad3579b7752fcc4254074dbea10
FIXTURE="${ROOT}/executions/borzoi-enformer-sequence-fixture-21064429"
FIXTURE_SHA256=d8a44d36c5b2b9997971e358e3553c10ce24c4badd41b28c91cfd48986123f0e
IMAGE="${RUNTIME}/tensorflow-tensorflow-2.17.0-gpu.sif"
TARGET="${ROOT}/executions/enformer-crested-restricted-forward-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
FIXTURE_ID=${FIXTURE_ID:?FIXTURE_ID must name one frozen Enformer fixture}
case "${FIXTURE_ID}" in
  trainfold2_ccre_032022|trainfold3_ccre_048015|trainfold4_ccre_064005) ;;
  *) printf 'Unsupported Enformer fixture: %s\n' "${FIXTURE_ID}" >&2; exit 1 ;;
esac

cleanup() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Enformer restricted forward: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "${RUNTIME_SHA256} ${RUNTIME}/ARTIFACTS.json" \
  "${ADMISSION_SHA256} ${ADMISSION}/ARTIFACTS.json" \
  "${INSPECTION_SHA256} ${INSPECTION}/ARTIFACTS.json" \
  "${CROSSWALK_SHA256} ${CROSSWALK}/ARTIFACTS.json" \
  "${FIXTURE_SHA256} ${FIXTURE}/ARTIFACTS.json"; do
  read -r digest path <<< "${pair}"
  printf '%s  %s\n' "${digest}" "${path}" | sha256sum --check --strict
done
install -d -m 0750 "${STAGE}/scratch" "${STAGE}/validation"

module load singularity/4.2.2
container=(
  singularity exec --nv --cleanenv
  --env PYTHONNOUSERSITE=1
  --env PYTHONDONTWRITEBYTECODE=1
  --env KERAS_BACKEND=tensorflow
  --env TF_DETERMINISTIC_OPS=1
  --env TF_CPP_MIN_LOG_LEVEL=1
  --bind "${ROOT}:${ROOT}"
  "${IMAGE}"
)
"${container[@]}" python "${ROOT}/tests/unit/enformer_crested_test_probe_checkpoint.py" \
  > "${STAGE}/validation/unittest.stdout.txt" \
  2> "${STAGE}/validation/unittest.stderr.txt"
"${container[@]}" python "${ROOT}/scripts/enformer_crested_probe_checkpoint.py" \
  --model-archive "${ADMISSION}/downloads/enformer_human.tar.gz" \
  --source-archive "${ADMISSION}/downloads/crested-source.tar.gz" \
  --fixture-fasta "${FIXTURE}/fixture/enformer.reference.fa.gz" \
  --fixture-manifest "${FIXTURE}/fixture/sequence_manifest.tsv" \
  --fixture-id "${FIXTURE_ID}" \
  --scratch "${STAGE}/scratch" \
  --output "${STAGE}/probe" \
  > "${STAGE}/validation/probe.stdout.txt"
rmdir "${STAGE}/scratch"
"${container[@]}" python --version > "${STAGE}/validation/python.version.txt"
"${container[@]}" python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${container[@]}" nvidia-smi -q > "${STAGE}/validation/nvidia-smi.txt"
module list 2> "${STAGE}/validation/module-list.txt"
sha256sum \
  "${ROOT}/scripts/enformer_crested_probe_checkpoint.py" \
  "${ROOT}/tests/unit/enformer_crested_test_probe_checkpoint.py" \
  "${ROOT}/slurm/enformer_crested_probe_checkpoint_l40s.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "probe/receipt.json").read_text())
if (
    receipt.get("status") != "pass"
    or receipt.get("native_sonnet_parity_established") is not False
    or receipt.get("outcomes_read") is not False
    or receipt.get("open_champion_eligible") is not False
):
    raise SystemExit("Enformer restricted checkpoint receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "restricted_sequence_model_checkpoint_forward",
        "model_id": receipt["model_id"],
        "checkpoint_member_sha256": receipt["checkpoint_member_sha256"],
        "fixture_id": receipt["fixture_id"],
        "parameter_count": receipt["parameter_count"],
        "deterministic_repeat_bit_identical": True,
        "native_sonnet_parity_established": False,
        "training_or_adaptation_executed": False,
        "outcomes_read": False,
        "open_champion_eligible": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Enformer CREsted restricted forward: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
