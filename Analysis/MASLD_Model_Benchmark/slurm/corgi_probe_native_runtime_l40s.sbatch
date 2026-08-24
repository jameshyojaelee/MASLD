#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-l40s-r8
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONHASHSEED=20260824

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-recovery-21064690"
RUNTIME_ARTIFACTS_SHA256=01de778c17c67d92de9a2e242350e17443dbfe76b7c5a5f8bc447eae5408cd3c
IMAGE="${RUNTIME}/nvidia-pytorch-24.01-py3.sif"
IMAGE_SHA256=eace72d1b831a4834d692134749bb82093c29930fb3a0258ae667b8acb840959
ACQUISITION="${ROOT}/executions/upstream-acquisitions/corgi-d84c000-zenodo18630048-21064463"
ACQUISITION_ARTIFACTS_SHA256=73aad9ab9ad3e080a2ccaa0cee1a1d727d126a454f0e69f74f343a900b168c0f
SCHEMA="${ROOT}/executions/corgi-official-checkpoint-schema-21064532"
SCHEMA_ARTIFACTS_SHA256=312bb3c90ae6da3eff20f43d920b03fbf6a5b05b2b480bbc98d24aaad71d868e
OVERLAY="${ROOT}/executions/environments/corgi-rotary-emb-0.1-flashattn2.0.4-build-r1-21064845"
OVERLAY_ARTIFACTS_SHA256=c60c4ebf0b32df7d158a3a1d781f5d45fdf4187ab0fba9f4bea35605a92e4c3a
ROTARY_PROBE="${ROOT}/executions/corgi-rotary-emb-l40s-21064869"
ROTARY_PROBE_ARTIFACTS_SHA256=4ac0d03f52fa26583006c4067ccf35ac4946a491295c721f5b6a7f4cbc90a3f2
SOURCE="${ACQUISITION}/source/corgi-d84c00082da3d04eb05169e290e53f1c7da2b9c5"
REGULAR="${ACQUISITION}/downloads/corgi_model.pt"
PLUS="${ACQUISITION}/downloads/corgiplus_model.pt"
PROBE_SHA256=f21d31c0c9888e3a80f730ff57ba4ccaf6b57cad887d542d16a442a241f34be6
TARGET="${ROOT}/executions/corgi-native-runtime-l40s-r8-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Corgi L40S probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" \
  "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${IMAGE_SHA256}" "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' "${ACQUISITION_ARTIFACTS_SHA256}" \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${SCHEMA_ARTIFACTS_SHA256}" \
  "${SCHEMA}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${OVERLAY_ARTIFACTS_SHA256}" \
  "${OVERLAY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${ROTARY_PROBE_ARTIFACTS_SHA256}" \
  "${ROTARY_PROBE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${PROBE_SHA256}" \
  "${ROOT}/scripts/corgi_probe_native_runtime.py" | sha256sum --check --strict
printf '%s  %s\n' 53f4891df15f2dfe204a6fbe133ae21e7dec983441ab5400c96360ff6965c665 \
  "${SOURCE}/corgi/modules.py" | sha256sum --check --strict
printf '%s  %s\n' 767067e4659ca196f307d2133104365fb7e6ea9d2f151dd61d3c2a8a7e0737e6 \
  "${SOURCE}/corgi/model.py" | sha256sum --check --strict
printf '%s  %s\n' a6ae55409602e0de35a8d2cf757eb5ebac93b602e31e00f0aa36a4fbd63ac4b1 \
  "${SOURCE}/corgi/config.py" | sha256sum --check --strict
printf '%s  %s\n' cd51539f01de1e66a3a4f9b89e772bdeae07df2f0178d5b692f2368db2f4d2a4 \
  "${REGULAR}" | sha256sum --check --strict
printf '%s  %s\n' 01855614bfeffa72d5bcaaa28729c022ccfbc8f5fdcc0add8e04ba4533e94ccd \
  "${PLUS}" | sha256sum --check --strict

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${RUNTIME}" "${ACQUISITION}" "${SCHEMA}" "${OVERLAY}" "${ROTARY_PROBE}" \
  "${STAGE}/validation/input-verification.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

roots = [Path(value).resolve(strict=True) for value in sys.argv[1:-1]]
if any("outcome" in part.lower() for root in roots for part in root.parts):
    raise SystemExit("evaluator outcome artifact entered the Corgi runtime probe")
payload = {root.name: verify_frozen_tree(root)["metadata"] for root in roots}
Path(sys.argv[-1]).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY

module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total \
  --format=csv,noheader > "${STAGE}/validation/nvidia-smi.csv"
singularity --version > "${STAGE}/validation/singularity.version.txt"
singularity exec --nv --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env PYTHONHASHSEED=20260824 \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  --env "PYTHONPATH=${OVERLAY}/site-packages" \
  "${IMAGE}" python "${ROOT}/scripts/corgi_probe_native_runtime.py" \
    --source "${SOURCE}" \
    --regular "${REGULAR}" \
    --plus "${PLUS}" \
    --output "${STAGE}/validation/probe.json" \
    > "${STAGE}/validation/probe.stdout.txt" \
    2> "${STAGE}/validation/probe.stderr.txt"

module purge
module load python/3.11.5-GCCcore-13.2.0
python - "${STAGE}/validation/probe.json" \
  "${STAGE}/validation/summary.json" <<'PY'
from pathlib import Path
import json
import sys

probe = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
if probe.get("status") != "pass" or probe.get("outcomes_exposed") is not False:
    raise SystemExit("Corgi L40S probe did not pass its firewall")
if probe.get("gpu_name") != "NVIDIA L40S" or probe.get("gpu_capability") != [8, 9]:
    raise SystemExit("Corgi L40S device contract differs")
for model_id in ("regular", "plus"):
    record = probe[model_id]
    if (
        record.get("strict_restore") is not True
        or record.get("finite_forward") is not True
        or record.get("finite_backward") is not True
        or record.get("full_524288_bp_forward_performed") is not True
        or record.get("full_524288_bp_output_shape") != [1, 22, 6144]
        or record.get("native_numeric_parity_performed") is not False
    ):
        raise SystemExit(f"Corgi runtime contract differs: {model_id}")
if probe["plus"].get("released_input_bundle_complete") is not False:
    raise SystemExit("Corgi+ incomplete released input bundle was obscured")
summary = {
    "schema_version": "masld-bench-corgi-native-runtime-l40s-v1",
    "campaign_revision": 8,
    "supersedes_failed_jobs": [21064720, 21064750, 21064883, 21064908, 21064928, 21064963, 21064998],
    "status": "pass",
    "regular_strict_restore": True,
    "plus_strict_restore": True,
    "regular_full_window_forward": True,
    "plus_full_window_forward": True,
    "shortened_window_backward": True,
    "native_numeric_parity_performed": False,
    "corgi_plus_released_input_bundle_complete": False,
    "evaluator_outcomes_exposed": False,
    "champion_eligibility_changed": False,
}
Path(sys.argv[2]).write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(summary, sort_keys=True))
PY

sha256sum \
  "${ROOT}/scripts/corgi_probe_native_runtime.py" \
  "${ROOT}/slurm/corgi_probe_native_runtime_l40s.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" \
  "${RUNTIME_ARTIFACTS_SHA256}" "${ACQUISITION_ARTIFACTS_SHA256}" \
  "${SCHEMA_ARTIFACTS_SHA256}" "${OVERLAY_ARTIFACTS_SHA256}" \
  "${ROTARY_PROBE_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

freeze_tree(
    Path(sys.argv[1]),
    {
        "artifact_class": "corgi_native_runtime_l40s_probe",
        "campaign_revision": 8,
        "supersedes_failed_jobs": [21064720, 21064750, 21064883, 21064908, 21064928, 21064963, 21064998],
        "runtime_artifacts_sha256": sys.argv[2],
        "acquisition_artifacts_sha256": sys.argv[3],
        "schema_artifacts_sha256": sys.argv[4],
        "rotary_overlay_artifacts_sha256": sys.argv[5],
        "rotary_l40s_probe_artifacts_sha256": sys.argv[6],
        "regular_strict_restore": True,
        "plus_strict_restore": True,
        "full_524288_bp_forward": True,
        "shortened_window_backward": True,
        "native_numeric_parity_performed": False,
        "corgi_plus_released_input_bundle_complete": False,
        "evaluator_outcomes_exposed": False,
        "champion_eligibility_changed": False,
        "status": "pass",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Corgi native L40S runtime probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
