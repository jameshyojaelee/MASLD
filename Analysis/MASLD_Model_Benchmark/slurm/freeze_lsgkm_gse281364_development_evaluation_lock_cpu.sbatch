#!/usr/bin/env bash
#SBATCH --job-name=model-check-625
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
CONFIG="${ROOT}/config/lsgkm_gse281364_dinucleotide_development_evaluation.json"
SCRIPT="${ROOT}/scripts/freeze_lsgkm_gse281364_development_evaluation_lock.py"
TARGET="${ROOT}/executions/model-check-625-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite evaluation lock: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' e9678f9c4c2919a65d67167f441b3938c5ad972b29b870269989969cfd8c0343 "${CONFIG}" | sha256sum --check --strict
printf '%s  %s\n' 144efe5e0c73db2108dd1f47c593fa7f17bf98c17329795449f0cc157eedf9d3 "${SCRIPT}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}" python "${SCRIPT}" \
  --project-root "${ROOT}" --config "${CONFIG}" \
  --output "${STAGE}/lock" \
  > "${STAGE}/lock.stdout.json" 2> "${STAGE}/lock.stderr.txt"
python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
sha256sum "${SCRIPT}" "${CONFIG}" \
  "${ROOT}/slurm/freeze_lsgkm_gse281364_development_evaluation_lock_cpu.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "lock/selection_lock.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "locked_before_evaluator_control_and_outcome_access"
    or len(receipt["candidates"]) != 3
    or receipt["prediction_values_opened"]
    or receipt["control_prediction_values_opened"]
    or receipt["outcomes_opened"]
    or receipt["metrics_calculated"]
    or receipt["sealed_assets_opened"]
):
    raise SystemExit("evaluation selection lock differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_development_evaluation_selection_lock",
    **receipt,
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'LS-GKM development evaluation lock: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
