#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-tile-preflight
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
WRAPPER="${ROOT}/slurm/run_corgi_hepatocyte_tile_smoke_bundle.sbatch"
TARGET="${ROOT}/executions/corgi-hepatocyte-tile-bundle-v3-preflight-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi bundle preflight artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' \
  3ee487623986aaffd999b85edad848a141c4a6f0151a5f2e2408bb0b9724da8f \
  "${WRAPPER}" | sha256sum --check --strict
bash -n "${WRAPPER}"
PREFLIGHT_ONLY=1 bash "${WRAPPER}" \
  > "${STAGE}/validation/preflight.stdout.txt" \
  2> "${STAGE}/validation/preflight.stderr.txt"
module purge
module load singularity/4.2.2
IMAGE="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-recovery-21064690/nvidia-pytorch-24.01-py3.sif"
ROOT_PROBE="${ROOT}/scripts/.corgi-read-only-probe-${SLURM_JOB_ID}"
STAGE_PROBE="${STAGE}/validation/nested-bind-write-probe.txt"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --bind "${STAGE}:${STAGE}:rw" \
  "${IMAGE}" python - "${ROOT_PROBE}" "${STAGE_PROBE}" <<'PY'
from pathlib import Path
import sys

read_only_probe = Path(sys.argv[1])
write_probe = Path(sys.argv[2])
try:
    read_only_probe.write_text("must fail\n")
except OSError:
    pass
else:
    raise SystemExit("benchmark root remained writable inside model container")
write_probe.write_text("nested staging bind is writable\n")
PY
if [[ -e "${ROOT_PROBE}" || ! -s "${STAGE_PROBE}" ]]; then
  printf 'Corgi model-container bind isolation differs.\n' >&2
  exit 1
fi
sha256sum \
  "${WRAPPER}" \
  "${ROOT}/slurm/test_corgi_hepatocyte_tile_bundle_preflight_cpu.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
stdout = (root / "validation/preflight.stdout.txt").read_text()
if "Corgi hepatocyte tile bundle preflight passed" not in stdout:
    raise SystemExit("Corgi tile bundle v3 preflight did not pass")
freeze_tree(
    root,
    {
        "artifact_class": "Corgi_outcome_aligned_tile_smoke_bundle_v3_preflight",
        "gpu_used": False,
        "outcomes_read": False,
        "wrapper_sha256": "3ee487623986aaffd999b85edad848a141c4a6f0151a5f2e2408bb0b9724da8f",
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi tile bundle v3 CPU preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
