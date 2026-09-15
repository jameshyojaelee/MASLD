#!/usr/bin/env bash
#SBATCH --job-name=model-check-624
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
CONFIG="${ROOT}/config/lsgkm_gse281364_dinucleotide_aggregation.json"
SCRIPT="${ROOT}/scripts/audit_lsgkm_gse281364_dinucleotide_aggregation.py"
AGGREGATION="${ROOT}/executions/model-data-623-21100461/ARTIFACTS.json"
AGGREGATION_SHA256=2d74167f7670bbcad7ec6e178f228b9f7efc80b308f827041d4df7df18357d1a
TARGET="${ROOT}/executions/model-check-624-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite LS-GKM aggregation audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' 5b801774fdb6ad6b652506797265a2bae352acbc7acd2eae8a6450e884ea482c "${SCRIPT}" | sha256sum --check --strict
printf '%s  %s\n' 8c1f2f2b07bff86d9fd5db2656819cccf19cad7ce6ad4ebb01ad8d145d6795d9 "${CONFIG}" | sha256sum --check --strict
printf '%s  %s\n' "${AGGREGATION_SHA256}" "${AGGREGATION}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}" python "${SCRIPT}" \
  --project-root "${ROOT}" --config "${CONFIG}" \
  --aggregation-artifacts "${AGGREGATION}" \
  --aggregation-artifacts-sha256 "${AGGREGATION_SHA256}" \
  --output "${STAGE}/audit" \
  > "${STAGE}/audit.stdout.json" 2> "${STAGE}/audit.stderr.txt"
python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
sha256sum "${SCRIPT}" "${CONFIG}" "${AGGREGATION}" \
  "${ROOT}/slurm/audit_lsgkm_gse281364_dinucleotide_aggregation_cpu.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "audit/aggregation_audit.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_independent_raw_aggregation_rederivation"
    or receipt["fit_count"] != 25
    or receipt["seed_by_element_scores"] != 5165
    or receipt["prediction_rows"] != 20660
    or receipt["context_copy_groups"] != 10330
    or receipt["control_prediction_member_opened"]
    or receipt["outcomes_read"]
    or receipt["metrics_calculated"]
    or receipt["sealed_assets_read"]
):
    raise SystemExit("LS-GKM aggregation audit receipt differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_exact_dinucleotide_aggregation_audit",
    **receipt,
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'LS-GKM aggregation audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
