#!/usr/bin/env bash
#SBATCH --job-name=model-work-151
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export MASLD_REQUIRE_CUDA=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/scvi-tools-1.4.2-inductive-portable-21079606
RUNTIME_ARTIFACTS_SHA256=d885ada24fdad9f6e3d57a92ff870ca8ff05f07701ac6b0762f55a8cb19db61a
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_ARTIFACTS_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_ARTIFACTS_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
EXPECTED_WHEEL_SHA256=c453b75b0fa0d222bf0eb3a531607d17f5ea0dbe5b848ab12fd8cc24599884e5
FIT_SCRIPT_SHA256=a3ecfba0cdf7b5ab912398ac3699c8f78f26906628e50fe7608376722d416e91
AGGREGATE_SCRIPT_SHA256=a465206b45e3a8f3230d6a953b1f10c98fed48c1b1899f15540d50b2038a0524
ADAPTER_SHA256=2e191ea97edd7af6a6e7ec75e6a312c168ed4fac7ef5ff7731981c2b9ea5c865

TARGET=${ROOT}/executions/model-work-151-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite scVI/scANVI seed bundle.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SPLIT_ARTIFACTS_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FIT_SCRIPT_SHA256}" "${ROOT}/scripts/fit_predict_scvi_scanvi_study_50000.py" \
  | sha256sum --check --strict
printf '%s  %s\n' "${AGGREGATE_SCRIPT_SHA256}" "${ROOT}/scripts/aggregate_scvi_scanvi_study_50000.py" \
  | sha256sum --check --strict
printf '%s  %s\n' "${ADAPTER_SHA256}" "${ROOT}/src/masld_bench/adapters/scvi_scanvi_inductive.py" \
  | sha256sum --check --strict

module purge
module load python/3.12.3-GCCcore-13.3.0
cd "${ROOT}"
export PYTHONPATH=${RUNTIME}/site-packages:${ROOT}:${ROOT}/src
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
python - "${EXPECTED_WHEEL_SHA256}" "${STAGE}" <<'PY' \
  > "${STAGE}/validation/gpu-runtime-probe.json"
import json
from importlib.metadata import version
import os
from pathlib import Path
import sys

import scvi
import torch

# Torch layer guard.  These runtimes are overlays stacked by PYTHONPATH, so a
# second complete torch install can sit one path element away and be picked up
# silently.  Prove the intended layer loaded, that it is the pinned build, and
# that the build's compiled architecture list actually covers this GPU, before
# any model code runs.  The result is written to validation/layer_guard.json
# whether it passes or fails.
EXPECTED_LAYER = "scvi-tools-1.4.2-inductive-portable-21079606"
EXPECTED_TORCH = "2.13.0+cu130"
EXPECTED_ARCH_LIST = ["sm_75", "sm_80", "sm_86", "sm_90", "sm_100", "sm_120"]
torch_file = os.path.realpath(torch.__file__)
arch_list = list(torch.cuda.get_arch_list())
guard = {
    "arch_list": arch_list,
    "cell_id": "scvi-scanvi-50k-study-held-3seed-v1",
    "cuda_available": torch.cuda.is_available(),
    "cuda_runtime": torch.version.cuda,
    "expected_layer": EXPECTED_LAYER,
    "expected_torch_version": EXPECTED_TORCH,
    "torch_file": torch_file,
    "torch_version": torch.__version__,
}
if guard["cuda_available"]:
    guard_capability = torch.cuda.get_device_capability(torch.device("cuda:0"))
    guard["device_capability"] = list(guard_capability)
    guard["device_count"] = torch.cuda.device_count()
    guard["device_name"] = torch.cuda.get_device_name(torch.device("cuda:0"))
    guard["device_sm"] = f"sm_{guard_capability[0]}{guard_capability[1]}"
    guard["device_sm_in_arch_list"] = guard["device_sm"] in arch_list
guard["layer_guard_passed"] = (
    EXPECTED_LAYER in torch_file
    and torch.__version__ == EXPECTED_TORCH
    and arch_list == EXPECTED_ARCH_LIST
    and guard["cuda_available"]
    and guard.get("device_count") == 1
    and guard.get("device_sm_in_arch_list") is True
)
Path(sys.argv[2], "validation/layer_guard.json").write_text(
    json.dumps(guard, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
if not guard["layer_guard_passed"]:
    raise SystemExit("torch layer guard failed: " + json.dumps(guard, sort_keys=True))

if (
    scvi.__version__ != "1.4.2"
    or version("scvi-tools") != "1.4.2"
    or not torch.cuda.is_available()
    or torch.cuda.device_count() != 1
):
    raise SystemExit("scvi-tools or one-GPU runtime differs")
device = torch.device("cuda:0")
capability = torch.cuda.get_device_capability(device)
left = torch.arange(256 * 64, dtype=torch.float32, device=device).reshape(256, 64)
right = torch.arange(64 * 32, dtype=torch.float32, device=device).reshape(64, 32)
result = left @ right
torch.cuda.synchronize()
if capability != (12, 0) or not torch.isfinite(result).all():
    raise SystemExit("B6K capability or CUDA kernel differs")
print(json.dumps({
    "status": "pass",
    "scvi_tools": scvi.__version__,
    "scvi_tools_wheel_sha256": sys.argv[1],
    "torch": torch.__version__,
    "torch_cuda_build": torch.version.cuda,
    "gpu_name": torch.cuda.get_device_name(device),
    "compute_capability": list(capability),
    "device_sm": guard["device_sm"],
    "forward_kernel": "pass",
    "layer_guard_passed": guard["layer_guard_passed"],
    "project_data_read": False,
    "sealed_outcomes_read": False,
}, sort_keys=True, separators=(",", ":")))
PY
nvidia-smi -q > "${STAGE}/validation/nvidia-smi-q.txt"
python -m unittest \
  tests.unit.test_scvi_scanvi_inductive \
  tests.unit.test_fit_predict_scvi_scanvi_study_50000 \
  tests.unit.test_aggregate_scvi_scanvi_study_50000 \
  > "${STAGE}/validation/unit-tests.txt" \
  2>&1
for SEED in 20260824 20260825 20260826; do
  python scripts/fit_predict_scvi_scanvi_study_50000.py \
    --source "${SOURCE}" \
    --split "${SPLIT}" \
    --output "${STAGE}/seed-${SEED}" \
    --expected-source-artifacts-sha256 "${SOURCE_ARTIFACTS_SHA256}" \
    --expected-split-artifacts-sha256 "${SPLIT_ARTIFACTS_SHA256}" \
    --seed "${SEED}" \
    > "${STAGE}/validation/seed-${SEED}.stdout.json" \
    2> "${STAGE}/validation/seed-${SEED}.stderr.txt"
done
python scripts/aggregate_scvi_scanvi_study_50000.py \
  --seed-root "${STAGE}/seed-20260824" \
  --seed-root "${STAGE}/seed-20260825" \
  --seed-root "${STAGE}/seed-20260826" \
  --output "${STAGE}/ensemble" \
  > "${STAGE}/validation/aggregate.stdout.json" \
  2> "${STAGE}/validation/aggregate.stderr.txt"
cp "${RUNTIME}/validation/pip-freeze.txt" "${STAGE}/validation/runtime-pip-freeze.txt"
sha256sum \
  scripts/fit_predict_scvi_scanvi_study_50000.py \
  scripts/aggregate_scvi_scanvi_study_50000.py \
  scripts/score_scvi_scanvi_study_50000.py \
  src/masld_bench/adapters/scvi_scanvi_inductive.py \
  tests/unit/test_scvi_scanvi_inductive.py \
  tests/unit/test_fit_predict_scvi_scanvi_study_50000.py \
  tests/unit/test_aggregate_scvi_scanvi_study_50000.py \
  slurm/run_scvi_scanvi_study_50000_seed_bundle.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${STAGE}" "${RUNTIME_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
for seed in (20260824, 20260825, 20260826):
    seed_root = root / f"seed-{seed}"
    receipt = json.loads(
        (seed_root / "prediction_receipt.json").read_text(encoding="utf-8")
    )
    if (
        receipt.get("status") != "pass_development_predictions"
        or receipt.get("seed") != seed
        or receipt.get("rows") != 50_000
        or receipt.get("donors") != 102
        or receipt.get("studies") != 7
        or receipt.get("metrics_calculated") is not False
        or receipt.get("common_head_fit") is not False
        or receipt.get("held_query_training") is not False
        or receipt.get("held_query_adaptation") is not False
        or receipt.get("sealed_outcomes_read") is not False
    ):
        raise SystemExit("scVI/scANVI prediction receipt differs")
    freeze_tree(
        seed_root,
        {
            "artifact_class": "scvi_scanvi_task_native_seed_predictions",
            "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
            "split_id": "resource_atlas_study_outer_5fold_v1",
            "seed": seed,
            "rows": 50_000,
            "donors": 102,
            "studies": 7,
            "models": ["scvi_baseline", "scanvi_baseline"],
            "query_training_executed": False,
            "query_adaptation_executed": False,
            "metrics_calculated": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )
    verify_frozen_tree(seed_root)
ensemble = json.loads(
    (root / "ensemble/ensemble_receipt.json").read_text(encoding="utf-8")
)
if (
    ensemble.get("status") != "pass_development_predictions"
    or ensemble.get("rows") != 50_000
    or ensemble.get("seeds") != [20260824, 20260825, 20260826]
    or ensemble.get("metrics_calculated") is not False
    or ensemble.get("held_query_training") is not False
    or ensemble.get("held_query_adaptation") is not False
    or ensemble.get("sealed_outcomes_read") is not False
):
    raise SystemExit("scVI/scANVI ensemble receipt differs")
freeze_tree(
    root / "ensemble",
    {
        "artifact_class": "scvi_scanvi_task_native_ensemble_predictions",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "split_id": "resource_atlas_study_outer_5fold_v1",
        "seeds": [20260824, 20260825, 20260826],
        "rows": 50_000,
        "models": ["scvi_baseline", "scanvi_baseline"],
        "metrics_calculated": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
verify_frozen_tree(root / "ensemble")
freeze_tree(
    root,
    {
        "artifact_class": "scvi_scanvi_task_native_three_seed_bundle",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "split_id": "resource_atlas_study_outer_5fold_v1",
        "seeds": [20260824, 20260825, 20260826],
        "runtime_artifacts_sha256": sys.argv[2],
        "gpu_runtime_probe": "passed",
        "metrics_calculated": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum \
  "${TARGET}/ARTIFACTS.json" \
  "${TARGET}/ensemble/ARTIFACTS.json" \
  "${TARGET}/seed-20260824/ARTIFACTS.json" \
  "${TARGET}/seed-20260825/ARTIFACTS.json" \
  "${TARGET}/seed-20260826/ARTIFACTS.json"
