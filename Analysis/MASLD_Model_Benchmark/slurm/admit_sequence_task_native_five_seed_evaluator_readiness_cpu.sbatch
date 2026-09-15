#!/usr/bin/env bash
#SBATCH --job-name=model-check-293
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${SLURM_JOB_ID:?SLURM_JOB_ID is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PLAN="${ROOT}/executions/sequence-five-seed-rectangle-admission-21098556"
CORRECTION="${ROOT}/executions/sequence-five-seed-execution-correction-21098951"
QUEUE="${ROOT}/executions/sequence-five-seed-queue-admission-v3-21098963"
CONTRACT="${ROOT}/config/evaluation/sequence_task_native_five_seed_rectangular_evaluator_20260825.json"
TARGET="${ROOT}/executions/sequence-five-seed-evaluator-readiness-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite five-seed evaluator readiness artifact.\n' >&2
  exit 2
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' 78d785157d8b11b38d18f99362b98bdec2a894e6455f42ac54cf681958da0b28 \
  "${PLAN}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 5f3c122972f86a16e615c52c98e1f20af15fc6dcbf2c83b4b0f53730b0a68176 \
  "${CORRECTION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' ceca0985719a13b7588dbc5156401f30b611906338aecc0bfaa07ce0eba5a646 \
  "${QUEUE}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}/src:${ROOT}/scripts" python -m unittest -v \
  tests.unit.test_sequence_task_native_five_seed_evaluator_readiness \
  tests.unit.test_stats.TwoWayBootstrapTests \
  > "${STAGE}/validation/unit.stdout.txt" \
  2> "${STAGE}/validation/unit.stderr.txt"
PYTHONPATH="${ROOT}/src" python \
  "${ROOT}/scripts/audit_sequence_task_native_five_seed_evaluator_readiness.py" \
  --root "${ROOT}" \
  --contract "${CONTRACT}" \
  --plan-root "${PLAN}" \
  --correction-root "${CORRECTION}" \
  --queue-root "${QUEUE}" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONTRACT}" \
  "${ROOT}/scripts/audit_sequence_task_native_five_seed_evaluator_readiness.py" \
  "${ROOT}/tests/unit/test_sequence_task_native_five_seed_evaluator_readiness.py" \
  "${ROOT}/src/masld_bench/evaluators/stats.py" \
  "${ROOT}/config/evaluation/inference.toml" \
  "${ROOT}/slurm/admit_sequence_task_native_five_seed_evaluator_readiness_cpu.sbatch" \
  > "${STAGE}/validation/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "audit/readiness_lock.json").read_text(encoding="utf-8"))
if (
    summary["prediction_or_metric_values_read"] is not False
    or summary["observed_outcome_files_read"] is not False
    or summary["partial_cross_model_ranking_authorized"] is not False
    or summary["cross_model_ranking_performed"] is not False
    or summary["test_role_authorized"] is not False
    or summary["sealed_data_authorized"] is not False
):
    raise SystemExit("readiness firewall differs")
freeze_tree(root, {
    "artifact_class": "sequence_task_native_five_seed_evaluator_readiness",
    "status": "passed",
    "gate_status": summary["status"],
    "expected_fits": summary["expected_fits"],
    "expected_prediction_views": summary["expected_prediction_views"],
    "fit_disposition_counts": summary["fit_disposition_counts"],
    "prediction_view_disposition_counts": summary["prediction_view_disposition_counts"],
    "full_compatible_rectangle": summary["full_compatible_rectangle"],
    "outcome_access_authorized": summary["outcome_access_authorized"],
    "development_cross_model_ranking_authorized": summary["development_cross_model_ranking_authorized"],
    "partial_cross_model_ranking_authorized": False,
    "cross_model_ranking_performed": False,
    "prediction_or_metric_values_read": False,
    "observed_outcome_files_read": False,
    "test_role_authorized": False,
    "sealed_data_authorized": False,
    "uncertainty": "paired_two_way_donor_genomic_block_bootstrap_10000",
    "seeds_are_biological_replicates": False,
})
PY
module purge
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
