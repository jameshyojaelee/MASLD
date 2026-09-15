#!/usr/bin/env bash
#SBATCH --job-name=model-check-288
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${SLURM_JOB_ID:?SLURM_JOB_ID is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/sequence-five-seed-rectangle-admission-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite rectangle admission artifact.\n' >&2
  exit 2
fi
install -d -m 0750 "${STAGE}/validation"

printf '%s  %s\n' 324de8ec6789c948380aaeb1fc7bb470083b21d4f77c3008c7b6870bd745390c \
  "${ROOT}/config/campaigns/sequence_task_native_five_seed_rectangle_20260825.json" | sha256sum --check --strict
printf '%s  %s\n' 1f1a09a372b0f4a2b023eb14912873106c1ddf48f2c7f093ce3fb69882ccdf8a \
  "${ROOT}/scripts/build_sequence_task_native_five_seed_rectangle.py" | sha256sum --check --strict
printf '%s  %s\n' 4da56214214fd2e876c9b4dd9d9d49746782f3a972fe62bf2bc7bd12208a0582 \
  "${ROOT}/scripts/run_sequence_task_native_five_seed_bundle.py" | sha256sum --check --strict
printf '%s  %s\n' 1bad33ca17916aa1fde2683044b85fb67fb8b107fc536368ecc1cee1b20aa78c \
  "${ROOT}/scripts/materialize_sequence_task_native_five_seed_queue.py" | sha256sum --check --strict
printf '%s  %s\n' d55fd67898661a8f652ee6ee74d5fcc07b01b4b03d96d08d1a79c9276f9295b2 \
  "${ROOT}/slurm/run_sequence_task_native_five_seed_bundle.sbatch" | sha256sum --check --strict
printf '%s  %s\n' ae55ff1d2068d43b9108149f3401fafa627f92878b3d824902b121b61e7cf3dc \
  "${ROOT}/tests/unit/test_sequence_task_native_five_seed_rectangle.py" | sha256sum --check --strict
printf '%s  %s\n' d22d3a840014921b7b8db1f801d37081f6e7092691228246c07cde7243c0386d \
  "${ROOT}/slurm/bpnet_train_predict_hepatocyte_donor0_genomic0.sbatch" | sha256sum --check --strict
printf '%s  %s\n' 81966ffb9d63f3f4ca902fd1cf9ee641ac30dbff2701fdafc1d24ae9402f7edd \
  "${ROOT}/slurm/train_chrombpnet_hepatocyte_donor0_genomic0_seed20260824.sbatch" | sha256sum --check --strict
printf '%s  %s\n' d0e256e09adb6416f7269493652fb5e1230e38a7e48c40a1ccba167fb9e103c3 \
  "${ROOT}/scripts/sequence_control_production_driver.sh" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python -m unittest -v \
  tests.unit.test_sequence_task_native_five_seed_rectangle \
  > "${STAGE}/validation/unit.stdout.txt" \
  2> "${STAGE}/validation/unit.stderr.txt"
PYTHONPATH="${ROOT}/src" python \
  "${ROOT}/scripts/build_sequence_task_native_five_seed_rectangle.py" \
  --root "${ROOT}" \
  --contract "${ROOT}/config/campaigns/sequence_task_native_five_seed_rectangle_20260825.json" \
  --output "${STAGE}/plan" \
  --frozen-root "${TARGET}" \
  --bundle-wrapper "${ROOT}/slurm/run_sequence_task_native_five_seed_bundle.sbatch"
PYTHONPATH="${ROOT}/src" python - "${STAGE}/plan" <<'PY'
from collections import Counter
import csv
import json
from pathlib import Path
import sys

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "summary.json").read_text(encoding="utf-8"))
estimate = json.loads((root / "resource_estimate.json").read_text(encoding="utf-8"))
with (root / "expected_fit_matrix.tsv").open(encoding="utf-8", newline="") as handle:
    matrix = list(csv.DictReader(handle, delimiter="\t"))
with (root / "bundle_tasks.tsv").open(encoding="utf-8", newline="") as handle:
    tasks = list(csv.DictReader(handle, delimiter="\t"))
with (root / "bundle_allocations.tsv").open(encoding="utf-8", newline="") as handle:
    bundles = list(csv.DictReader(handle, delimiter="\t"))
if (
    summary["status"] != "pass"
    or len(matrix) != 500
    or Counter(row["fit_state"] for row in matrix)["missing"] != len(tasks)
    or len(tasks) != summary["missing_fits"]
    or len(bundles) != summary["bundle_allocations"]
    or sum(int(row["logical_fits"]) for row in bundles) != len(tasks)
    or estimate["missing_fit_count"] != len(tasks)
    or summary["expected_prediction_views"] != 625
    or summary["prediction_or_metric_values_read"] is not False
    or summary["outcome_paths_available"] is not False
    or summary["manual_gpu_submission_allowed"] is not False
):
    raise SystemExit("rectangle admission invariant failed")
print(json.dumps({
    "status": "pass",
    "complete_reusable_fits": summary["complete_reusable_fits"],
    "missing_fits": summary["missing_fits"],
    "bundle_allocations": summary["bundle_allocations"],
    "estimated_storage_bytes": estimate["estimated_storage_bytes"],
    "median_estimated_gpu_hours": estimate["median_estimated_gpu_hours"],
    "p90_estimated_gpu_hours": estimate["p90_estimated_gpu_hours"],
}, sort_keys=True))
PY
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/config/campaigns/sequence_task_native_five_seed_rectangle_20260825.json" \
  "${ROOT}/scripts/build_sequence_task_native_five_seed_rectangle.py" \
  "${ROOT}/scripts/run_sequence_task_native_five_seed_bundle.py" \
  "${ROOT}/scripts/materialize_sequence_task_native_five_seed_queue.py" \
  "${ROOT}/slurm/run_sequence_task_native_five_seed_bundle.sbatch" \
  "${ROOT}/tests/unit/test_sequence_task_native_five_seed_rectangle.py" \
  "${ROOT}/slurm/bpnet_train_predict_hepatocyte_donor0_genomic0.sbatch" \
  "${ROOT}/slurm/train_chrombpnet_hepatocyte_donor0_genomic0_seed20260824.sbatch" \
  "${ROOT}/scripts/sequence_control_production_driver.sh" \
  "${ROOT}/slurm/bpnet_train_predict_five_seed_v2.sbatch" \
  "${ROOT}/slurm/chrombpnet_train_predict_five_seed_v2.sbatch" \
  "${ROOT}/scripts/sequence_control_production_driver_five_seed_v2.sh" \
  "${ROOT}/slurm/sequence_control_train_cnn_five_seed_v2.sbatch" \
  "${ROOT}/slurm/sequence_control_train_transformer_five_seed_v2.sbatch" \
  > "${STAGE}/validation/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "plan/summary.json").read_text(encoding="utf-8"))
estimate = json.loads((root / "plan/resource_estimate.json").read_text(encoding="utf-8"))
freeze_tree(root, {
    "artifact_class": "sequence_task_native_five_seed_rectangle_admission",
    "dataset_id": "gse296875",
    "task_id": "rna_conditioned_atac",
    "status": "passed",
    "expected_fits": 500,
    "expected_prediction_views": 625,
    "complete_reusable_fits": summary["complete_reusable_fits"],
    "missing_fits": summary["missing_fits"],
    "bundle_allocations": summary["bundle_allocations"],
    "estimated_storage_bytes": estimate["estimated_storage_bytes"],
    "median_estimated_gpu_hours": estimate["median_estimated_gpu_hours"],
    "p90_estimated_gpu_hours": estimate["p90_estimated_gpu_hours"],
    "full_rectangle_retained": True,
    "chrombpnet_views_share_one_fit": True,
    "prediction_or_metric_values_read": False,
    "outcome_paths_available": False,
    "sealed_paths_available": False,
    "cross_model_ranking_performed": False,
    "manual_gpu_submission_allowed": False,
    "qos": "nslab",
})
PY
module purge
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
