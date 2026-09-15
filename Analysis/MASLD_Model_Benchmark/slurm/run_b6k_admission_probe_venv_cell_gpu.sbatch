#!/usr/bin/env bash
#SBATCH --job-name=model-probe-010
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# One cell of the frozen B6K inclusion 2x2.  Submitted directly rather than through
# the GPU dependency dispatcher: the user authorized B6K in full on 2026-08-25, and
# the dispatcher hardcodes gpu:l40s:1 in four independent places so it physically
# cannot route a B6K job.  Job ids from this wrapper never enter the dispatcher's
# roster, which is built from its own state directory, so job 21100574 is unaffected.
#
# Synthetic tensors only.  Reads no biological data, no outcome, no held-back source,
# and calculates no benchmark metric.
#
# Required exports: B6K_CELL_ID, B6K_PYTHONPATH, B6K_TOLERANCE_SHA256

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TORCH_OVERLAY="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
COBOLT_CKPT="${ROOT}/executions/cobolt-runtime-probe-21060919/state_dict.pt"
COBOLT_SHA256=fa329ed1f4649c1cc375ca042211e09ee74402d4b1ba5df9aaca0a56edbadaeb
TOLERANCE="${ROOT}/config/b6k_numeric_tolerance_prespecification_20260825.json"

: "${B6K_CELL_ID:?B6K_CELL_ID is required}"
: "${B6K_PYTHON:?B6K_PYTHON is required}"
: "${B6K_TOLERANCE_SHA256:?B6K_TOLERANCE_SHA256 is required}"

TARGET="${ROOT}/executions/b6k-admission-cell-${B6K_CELL_ID}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite B6K cell artifact: %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}"

sha256sum --check --strict <<EOF
${COBOLT_SHA256}  ${COBOLT_CKPT}
${B6K_TOLERANCE_SHA256}  ${TOLERANCE}
4474a0d069d09a9a805e7d6e33936af95eeab9fb72d6e3c131e7e7733350842d  ${ROOT}/executions/environments/scprint2-native-torch2.8.0-v2-21073597/ARTIFACTS.json
EOF

cd "${ROOT}"
PY_BIN="${B6K_PYTHON}"
nvidia-smi > "${STAGE}/nvidia-smi.txt" 2>&1 || true

# Layer guard: prove the intended torch layer is the one that actually loaded.
# A complete torch 2.13.0+cu130 sits one PYTHONPATH element away in the overlay
# stack, held out only by element order, and nothing else asserts that order.
PYTHONPATH="${ROOT}" "${PY_BIN}" - "${STAGE}" "${B6K_CELL_ID}" <<'PY'
import json, sys, os
from pathlib import Path
stage, cell = Path(sys.argv[1]), sys.argv[2]
import torch
expected = "scprint2-native-torch2.8.0-v2-21073597"
resolved = os.path.realpath(torch.__file__)
guard = {"cell_id": cell, "expected_layer": expected, "torch_file": resolved,
         "torch_version": torch.__version__, "layer_guard_passed": expected in resolved,
         "arch_list": torch.cuda.get_arch_list()}
(stage / "layer_guard.json").write_text(json.dumps(guard, indent=2, sort_keys=True) + "\n")
print(json.dumps(guard, sort_keys=True))
if not guard["layer_guard_passed"]:
    raise SystemExit(f"layer guard failed: expected {expected} in {resolved}")
PY

PYTHONPATH="${ROOT}" "${PY_BIN}" -m unittest -v tests.unit.test_b6k_admission_probe \
  > "${STAGE}/unit_tests.stdout.txt" 2> "${STAGE}/unit_tests.stderr.txt"

set +e
PYTHONPATH="${ROOT}" "${PY_BIN}" -m scripts.b6k_admission_probe \
  --cell-id "${B6K_CELL_ID}" \
  --checkpoint "${COBOLT_CKPT}" \
  --checkpoint-sha256 "${COBOLT_SHA256}" \
  --output "${STAGE}/probe" \
  > "${STAGE}/probe.stdout.txt" 2> "${STAGE}/probe.stderr.txt"
PROBE_RC=$?
set -e
printf '%s\n' "${PROBE_RC}" > "${STAGE}/probe.returncode.txt"

python3 -m pip freeze --all > "${STAGE}/pip_freeze.txt" 2>/dev/null || true
sha256sum \
  "${ROOT}/slurm/run_b6k_admission_probe_venv_cell_gpu.sbatch" \
  "${ROOT}/scripts/b6k_admission_probe.py" \
  "${ROOT}/tests/unit/test_b6k_admission_probe.py" \
  "${TOLERANCE}" > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python - "${STAGE}" "${B6K_CELL_ID}" "${PROBE_RC}" <<'PY'
import json, sys
from pathlib import Path
from masld_bench.artifacts import freeze_tree
stage, cell, rc = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3])
probe = stage / "probe" / "probe.json"
status = "unknown"
if probe.exists():
    status = json.loads(probe.read_text()).get("status", "unknown")
freeze_tree(stage, {
    "artifact_class": "b6k_admission_probe_cell",
    "cell_id": cell,
    "probe_returncode": rc,
    "probe_status": status,
    "synthetic_inputs_only": True,
    "biological_data_read": False,
    "development_outcomes_read": False,
    "sealed_outcomes_read": False,
    "benchmark_metrics_calculated": False,
})
print(json.dumps({"cell_id": cell, "probe_status": status, "returncode": rc}, sort_keys=True))
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
