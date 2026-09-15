#!/usr/bin/env bash
#SBATCH --job-name=model-check-232
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
CONVERSION="${ROOT}/executions/environments/dna-lm-torch2.6-conversion-21064945"
CONFIG="${ROOT}/config/gse281364_task_native_sequence_control_production.toml"
TASK="${ROOT}/executions/model-check-227-21092007"
OUTCOMES="${ROOT}/executions/gse281364-replicate-outcomes-21066307"
TARGET="${ROOT}/executions/model-check-232-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed production preflight: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite production preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation/materialization"

printf '%s  %s\n' \
  622ca08fa8df12f91349431c40e82bde693c56423c755dea6f72dd76c5478a5f \
  "${TASK}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' \
  cb42bd4037127a4e4d6a406d6f7b00c706bc93b631cdf93083ce994809e791db \
  "${OUTCOMES}/ARTIFACTS.json" | sha256sum --check --strict
for binding in \
  "845820b7d7a510e6b959bb318914da161165fd71771f0a139b74c31da96bafd9 config/gse281364_task_native_sequence_control_production.toml" \
  "58b83a584fab05cded6f254817548125a3fb489abe230e783c7b1d88cbf7cb15 scripts/gse281364_task_native_sequence_controls.py" \
  "79add1142cdaa7a5b4b4849d33e26b78ea91b91b50a7f065dd2296cbe470ffa9 scripts/materialize_gse281364_sequence_control_training_targets.py" \
  "9b37c5050d36c394bcb3c916d4743e6e5ea7ffe3178ab211146f460e85e6eb07 scripts/train_predict_gse281364_task_native_sequence_controls.py" \
  "7c274121e1fb63d1db1c0bf53b2efa3ac19c4674728317345c3278f545a78535 scripts/aggregate_gse281364_task_native_sequence_control_predictions.py" \
  "93eb6889e025d28c9da6bd522f52b6c3c04570d154a2d1f21ee3c2bc3a6c6a49 tests/unit/test_gse281364_task_native_sequence_controls.py" \
  "f7e5efb483021a2f644c35ae54d6f758e6693979664f324d3eaf2c6da2e6e831 tests/unit/test_train_predict_gse281364_task_native_sequence_controls.py"; do
  read -r expected relative <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${TASK}" "${OUTCOMES}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
print("production input trees verified")
PY
module purge

module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${CONVERSION}/site-packages:${TORCH}/site-packages:${BASE}/site-packages:${ROOT}/src:${ROOT}"
python -m unittest -v \
  tests.unit.test_gse281364_task_native_sequence_controls \
  tests.unit.test_train_predict_gse281364_task_native_sequence_controls \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
python - "${CONFIG}" "${ROOT}/scripts/train_predict_gse281364_task_native_sequence_controls.py" \
  "${ROOT}/scripts/aggregate_gse281364_task_native_sequence_control_predictions.py" \
  > "${STAGE}/validation/static_contract.json" <<'PY'
from pathlib import Path
import json
import sys
import tomllib

config = tomllib.loads(Path(sys.argv[1]).read_text())
trainer = Path(sys.argv[2]).read_text().lower()
aggregate = Path(sys.argv[3]).read_text().lower()
if "tbd" in Path(sys.argv[1]).read_text().lower():
    raise SystemExit("production config retains TBD")
if any(token in trainer for token in (
    "replicate_outcomes", "spearman", "pearson", "r2_score",
    "roc_auc", "average_precision", "gse281364-replicate-outcomes",
)):
    raise SystemExit("training adapter source crosses outcome or metric firewall")
if any(token in aggregate for token in (
    "replicate_outcomes", "spearman", "pearson", "r2_score",
    "roc_auc", "average_precision",
)):
    raise SystemExit("aggregator source crosses outcome or metric firewall")
if (
    config["execution"]["fit_count"] != 50
    or config["execution"]["prediction_rows_total"] != 20660
    or config["output_contract"]["checkpoint_count"] != 50
    or config["output_contract"]["training_receipt_count"] != 50
    or config["output_contract"]["benchmark_metrics"] != "forbidden"
):
    raise SystemExit("production count or evaluator firewall differs")
print(json.dumps({
    "status": "pass",
    "fit_count": 50,
    "checkpoint_count": 50,
    "prediction_rows": 20660,
    "trainer_has_outcome_path": False,
    "trainer_has_benchmark_metrics": False,
    "aggregator_has_outcome_path": False,
    "aggregator_has_benchmark_metrics": False,
}, indent=2, sort_keys=True))
PY

PRIVATE="$(mktemp -d -p "${SLURM_TMPDIR:-/tmp}" sequence-target-preflight.XXXXXX)"
started="$(date +%s)"
for fold in fold-0 fold-1 fold-2 fold-3 fold-4; do
  python "${ROOT}/scripts/materialize_gse281364_sequence_control_training_targets.py" \
    --project-root "${ROOT}" \
    --config "${CONFIG}" \
    --outer-fold "${fold}" \
    --output "${PRIVATE}/${fold}" \
    > "${STAGE}/validation/materialization/${fold}.stdout.json" \
    2> "${STAGE}/validation/materialization/${fold}.stderr.txt"
  cp "${PRIVATE}/${fold}/receipt.json" \
    "${STAGE}/validation/materialization/${fold}.receipt.json"
done
python - "${STAGE}/validation/materialization" > "${STAGE}/validation/materialization_contract.json" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1])
held = {"fold-0": 208, "fold-1": 208, "fold-2": 207, "fold-3": 208, "fold-4": 202}
for fold, held_count in held.items():
    receipt = json.loads((root / f"{fold}.receipt.json").read_text())
    if (
        receipt["status"] != "pass_outer_training_targets_only"
        or receipt["outer_test_fold"] != fold
        or receipt["held_test_elements"] != held_count
        or receipt["training_elements"] != 1033 - held_count
        or receipt["training_target_rows"] != 2 * (1033 - held_count)
        or receipt["held_test_target_rows_emitted"] != 0
        or receipt["held_test_outcomes_visible_to_training_adapter"]
        or receipt["benchmark_metrics_calculated"]
        or receipt["raw_target_view_intended_for_frozen_output"]
    ):
        raise SystemExit("training-only materialization receipt differs")
print(json.dumps({
    "status": "pass",
    "outer_folds": 5,
    "held_test_target_rows_emitted": 0,
    "training_target_views_hash_bound": True,
}, indent=2, sort_keys=True))
PY
find "${PRIVATE}" -type f -delete
find "${PRIVATE}" -depth -type d -empty -delete
if [[ -e "${PRIVATE}" ]]; then
  printf 'Private target views remain after preflight cleanup\n' >&2
  exit 1
fi
ended="$(date +%s)"
python - "${started}" "${ended}" > "${STAGE}/validation/runtime_measurement.json" <<'PY'
import json
import sys
elapsed = int(sys.argv[2]) - int(sys.argv[1])
print(json.dumps({
    "status": "pass",
    "five_fold_real_target_materialization_elapsed_seconds": elapsed,
    "fixture_architecture_validation_elapsed_seconds": 1.208,
    "fixture_architecture_validation_source": "model-check-227-21092007",
    "production_device_recommendation": "one_L40S_due_to_100_training_stages_including_6_layer_230_token_self_attention",
}, indent=2, sort_keys=True))
PY
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/gse281364_task_native_sequence_controls.py" \
  "${ROOT}/scripts/materialize_gse281364_sequence_control_training_targets.py" \
  "${ROOT}/scripts/train_predict_gse281364_task_native_sequence_controls.py" \
  "${ROOT}/scripts/aggregate_gse281364_task_native_sequence_control_predictions.py" \
  "${ROOT}/tests/unit/test_gse281364_task_native_sequence_controls.py" \
  "${ROOT}/tests/unit/test_train_predict_gse281364_task_native_sequence_controls.py" \
  "${ROOT}/slurm/preflight_gse281364_sequence_control_production_cpu.sbatch" \
  > "${STAGE}/source.sha256"
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
static = json.loads((root / "validation/static_contract.json").read_text())
materialization = json.loads(
    (root / "validation/materialization_contract.json").read_text()
)
runtime = json.loads((root / "validation/runtime_measurement.json").read_text())
if (
    static["status"] != "pass"
    or materialization["status"] != "pass"
    or runtime["status"] != "pass"
    or any(path.name == "training_targets.tsv" for path in root.rglob("*"))
):
    raise SystemExit("production preflight differs")
freeze_tree(root, {
    "artifact_class": "gse281364_task_native_sequence_control_production_preflight",
    "dataset_id": "gse281364",
    "fit_count_planned": 50,
    "checkpoint_count_planned": 50,
    "prediction_rows_planned": 20660,
    "held_test_target_rows_emitted": 0,
    "raw_training_targets_retained": False,
    "benchmark_metrics_calculated": False,
    "gpu_jobs_submitted": 0,
    "status": "passed",
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Sequence-control production preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
