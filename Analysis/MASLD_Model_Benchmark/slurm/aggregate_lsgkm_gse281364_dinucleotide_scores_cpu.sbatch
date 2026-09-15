#!/usr/bin/env bash
#SBATCH --job-name=model-data-623
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
SCRIPT="${ROOT}/scripts/aggregate_lsgkm_gse281364_dinucleotide_scores.py"
TEST="${ROOT}/tests/unit/test_aggregate_lsgkm_gse281364_dinucleotide_scores.py"
INCIDENT="${ROOT}/config/artifacts/incidents/lsgkm_shared_control_prediction_exposure_20260825.json"
TARGET="${ROOT}/executions/model-data-623-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite LS-GKM raw aggregation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

printf '%s  %s\n' 8c1f2f2b07bff86d9fd5db2656819cccf19cad7ce6ad4ebb01ad8d145d6795d9 "${CONFIG}" | sha256sum --check --strict
printf '%s  %s\n' 3a0367a350893989228536478acd83dcd052063342b5f9fc5fb67195452efc89 "${SCRIPT}" | sha256sum --check --strict
printf '%s  %s\n' 7f9e0b9365409f16202b7de27ab892cd9bc0bedbd4f2f0183992b2543b047337 "${TEST}" | sha256sum --check --strict
printf '%s  %s\n' 8485d619d45f9a9d22297c179193364619a205f90f49efc0d4e6db26c494c314 "${INCIDENT}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}" python -m unittest -v \
  tests.unit.test_aggregate_lsgkm_gse281364_dinucleotide_scores \
  > "${STAGE}/tests.stdout.txt" 2> "${STAGE}/tests.stderr.txt"
PYTHONPATH="${ROOT}" python "${SCRIPT}" \
  --project-root "${ROOT}" --config "${CONFIG}" \
  --output "${STAGE}/aggregation" \
  > "${STAGE}/aggregation.stdout.json" 2> "${STAGE}/aggregation.stderr.txt"
python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
sha256sum "${SCRIPT}" "${TEST}" "${CONFIG}" "${INCIDENT}" \
  "${ROOT}/slurm/aggregate_lsgkm_gse281364_dinucleotide_scores_cpu.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "aggregation/aggregation_receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_outcome_blind_raw_score_aggregation"
    or receipt["fit_count"] != 25
    or receipt["seed_by_element_scores"] != 5165
    or receipt["prediction_rows_per_readout"] != 10330
    or receipt["total_prediction_rows"] != 20660
    or receipt["static_native_score_directly_stackable"]
    or receipt["shared_control_prediction_member_opened"]
    or receipt["benchmark_metrics_calculated"]
    or any(receipt[field] for field in ("outcomes_read", "reporter_counts_read", "control_prediction_values_read", "sealed_assets_read"))
):
    raise SystemExit("LS-GKM raw aggregation receipt differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_exact_dinucleotide_raw_score_aggregation",
    "design_id": receipt["design_id"],
    "fit_count": receipt["fit_count"],
    "prediction_rows_per_readout": receipt["prediction_rows_per_readout"],
    "total_prediction_rows": receipt["total_prediction_rows"],
    "raw_predictions_sha256": receipt["raw_predictions_sha256"],
    "static_native_score_directly_stackable": False,
    "shared_control_prediction_member_opened": False,
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
printf 'LS-GKM raw aggregation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
