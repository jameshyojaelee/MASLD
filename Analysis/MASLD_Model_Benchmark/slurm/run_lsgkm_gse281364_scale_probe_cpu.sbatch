#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-train-612
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG="${ROOT}/config/lsgkm_gse281364_dinucleotide_scale_probe.json"
SCRIPT="${ROOT}/scripts/run_lsgkm_gse281364_scale_probe.py"
ADAPTER="${ROOT}/scripts/lsgkm_safe_adapter.py"
TEST="${ROOT}/tests/unit/test_run_lsgkm_gse281364_scale_probe.py"
TARGET="${ROOT}/executions/model-cpu-train-612-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite LS-GKM scale probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

printf '%s  %s\n' 7061c9cfbef02a539c40222fad13eb20e781e6ee6711d9755ec05c3cab9bc224 "${CONFIG}" | sha256sum --check --strict
printf '%s  %s\n' 57f81394dc4d6b91cbf69809c47c185719cc24521c439680800477109411560f "${SCRIPT}" | sha256sum --check --strict
printf '%s  %s\n' bef248a9f9c084e503b24d4f1dfb91774287a84ffd5cf1afa2b1c1c97e48596f "${ADAPTER}" | sha256sum --check --strict
printf '%s  %s\n' 59f812c30dcfbc6af37f378513af43e1d42a2bb4eec38a6cce39fcddabf3ddeb "${TEST}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}" python -m unittest -v \
  tests.unit.test_run_lsgkm_gse281364_scale_probe \
  > "${STAGE}/tests.stdout.txt" 2> "${STAGE}/tests.stderr.txt"
PYTHONPATH="${ROOT}" python "${SCRIPT}" \
  --project-root "${ROOT}" --config "${CONFIG}" \
  --output "${STAGE}/probe" \
  > "${STAGE}/probe.stdout.json" 2> "${STAGE}/probe.stderr.txt"
python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
sha256sum "${SCRIPT}" "${ADAPTER}" "${TEST}" "${CONFIG}" \
  "${ROOT}/slurm/run_lsgkm_gse281364_scale_probe_cpu.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "probe/scale_probe_receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_one_fit_production_command_scale_probe"
    or receipt["positive_records"] != 10000
    or receipt["negative_records"] != 10000
    or receipt["held_elements_scored"] <= 0
    or receipt["required_canonical_11mers_scored"] <= 0
    or receipt["maximum_resident_set_kbytes"] <= 0
    or receipt["production_fits_executed"] != 1
    or receipt["benchmark_metrics_calculated"]
    or any(receipt[field] for field in ("outcomes_read", "reporter_counts_read", "control_prediction_values_read", "sealed_assets_read"))
):
    raise SystemExit("LS-GKM scale probe receipt differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_exact_dinucleotide_null_scale_probe",
    "design_id": receipt["design_id"],
    "split_id": receipt["active_fit"]["split_id"],
    "model_seed": receipt["active_fit"]["model_seed"],
    "held_elements_scored": receipt["held_elements_scored"],
    "maximum_resident_set_kbytes": receipt["maximum_resident_set_kbytes"],
    "total_native_wall_seconds": receipt["total_native_wall_seconds"],
    "production_fits_executed": 1,
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
printf 'LS-GKM scale probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
