#!/usr/bin/env bash
#SBATCH --job-name=b6kpreflight
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# CPU-only preflight for the B6K (sm_120) inclusion question.  It submits no GPU
# work and touches no dispatcher file.  It answers three things the GPU probes
# should not have to discover for themselves:
#
#   1. Which torch actually wins under the campaign's stacked PYTHONPATH, given
#      that the cu130 base layer also ships a torch.
#   2. Whether the CUDA 13 packages that appear in the cu124 overlay's
#      effective-pip-freeze are reachable at run time, or are a freeze-time
#      union artifact of the overlay stack.
#   3. Whether the candidate probe checkpoint restores under weights_only=True.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TORCH_OVERLAY="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
CU130="${ROOT}/executions/environments/scvi-tools-1.4.2-inductive-portable-21079606"
COBOLT="${ROOT}/executions/cobolt-runtime-probe-21060919"
RNA_ATAC="${ROOT}/executions/masld-rna-atac-torch-probe-21056987"

TARGET="${ROOT}/executions/b6k-preflight-runtime-resolution-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite preflight artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

# Read-only inputs, verified before anything is imported.
sha256sum --check --strict <<EOF
30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af  ${BASE}/ARTIFACTS.json
d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3  ${TORCH_OVERLAY}/ARTIFACTS.json
fa329ed1f4649c1cc375ca042211e09ee74402d4b1ba5df9aaca0a56edbadaeb  ${COBOLT}/state_dict.pt
0f7d8b232f2430abfbd707670733c0d627497f6061ab4e2f0f65002d361fccb6  ${RNA_ATAC}/state_dict.pt
EOF

module load python/3.12.3-GCCcore-13.3.0

# Cell A reproduces exactly what every current cu124 job puts on sys.path.
# Cell B is the candidate cu130 runtime for the B6K arm, used standalone.
run_cell() {
  local label="$1" pythonpath="$2"
  PYTHONPATH="${pythonpath}:${ROOT}" python "${ROOT}/scripts/inspect_torch_runtime_resolution.py" \
    --label "${label}" \
    --checkpoint "${COBOLT}/state_dict.pt" \
    --checkpoint "${RNA_ATAC}/state_dict.pt" \
    --output "${STAGE}/resolution-${label}.json"
  # The overlay layer ships no pip of its own; it borrows the base layer's.
  # A missing pip is worth recording but is not a reason to fail the preflight.
  PYTHONPATH="${pythonpath}:${ROOT}" python -m pip freeze --all \
    > "${STAGE}/pip-freeze-${label}.txt" \
    2> "${STAGE}/pip-freeze-${label}.stderr.txt" \
    || printf 'pip freeze unavailable for %s\n' "${label}" \
       >> "${STAGE}/pip-freeze-${label}.stderr.txt"
}

run_cell "cu124-overlay-stack" "${TORCH_OVERLAY}/site-packages:${BASE}/site-packages"
run_cell "cu130-standalone" "${CU130}/site-packages"

# The probe suite's own unit tests, run inside the frozen campaign runtime so a
# torch-version incompatibility in the probe code shows up here and not on a GPU.
for cell in "cu124-overlay-stack:${TORCH_OVERLAY}/site-packages:${BASE}/site-packages" \
            "cu130-standalone:${CU130}/site-packages"; do
  label="${cell%%:*}"
  path="${cell#*:}"
  PYTHONPATH="${path}:${ROOT}" python -m unittest -v \
    tests.unit.test_b6k_admission_probe \
    > "${STAGE}/unittest-${label}.stdout.txt" \
    2> "${STAGE}/unittest-${label}.stderr.txt"
done

sha256sum \
  "${ROOT}/scripts/b6k_admission_probe.py" \
  "${ROOT}/scripts/inspect_torch_runtime_resolution.py" \
  "${ROOT}/tests/unit/test_b6k_admission_probe.py" \
  "${ROOT}/slurm/preflight_b6k_runtime_resolution_cpu.sbatch" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    stage,
    {
        "artifact_class": "b6k_preflight_runtime_resolution",
        "gpu_submitted": False,
        "dispatcher_files_touched": False,
        "outcomes_read": False,
        "project_data_read": False,
        "metrics_calculated": False,
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'preflight artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
