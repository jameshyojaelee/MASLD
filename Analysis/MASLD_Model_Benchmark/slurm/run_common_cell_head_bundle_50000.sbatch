#!/usr/bin/env bash
#SBATCH --job-name=model-training
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONHASHSEED=20260824
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
: "${MASLD_GPU_DISPATCHER_AUTHORITY:?dispatcher authority is required}"
if [[ "${MASLD_GPU_DISPATCHER_AUTHORITY}" != sequence_gpu_dispatcher_v2 ]]; then
  printf 'Invalid GPU dispatcher authority.\n' >&2
  exit 1
fi
: "${EMBEDDINGS_ROOT:?frozen embedding root is required}"
: "${EMBEDDINGS_ARTIFACTS_SHA256:?frozen embedding ARTIFACTS hash is required}"
: "${EMBEDDINGS_RELATIVE:?relative embedding bundle is required}"
: "${HEAD_IDS:=linear two_layer_mlp}"
: "${SCREEN_SEEDS:=1103 1201 1301}"
: "${OUTER_FOLDS:=0 1 2 3 4}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
TARGET=${ROOT}/executions/model-training-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite common-head bundle.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/shards" "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SPLIT_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${EMBEDDINGS_ARTIFACTS_SHA256}" "${EMBEDDINGS_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" - <<'PY' > "${STAGE}/validation/device.json"
import json
import os
import torch

value = {
    "capability": list(torch.cuda.get_device_capability(0)) if torch.cuda.is_available() else None,
    "count": torch.cuda.device_count(),
    "cublas_workspace_config": os.environ.get("CUBLAS_WORKSPACE_CONFIG"),
    "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
    "python_hash_seed": os.environ.get("PYTHONHASHSEED"),
}
if (
    value["count"] != 1
    or "L40S" not in str(value["device"])
    or value["capability"] != [8, 9]
    or value["cublas_workspace_config"] != ":4096:8"
    or value["python_hash_seed"] != "20260824"
):
    raise SystemExit("common-head GPU runtime differs")
print(json.dumps(value, sort_keys=True))
PY

read -r -a HEAD_ARRAY <<< "${HEAD_IDS}"
read -r -a SEED_ARRAY <<< "${SCREEN_SEEDS}"
read -r -a FOLD_ARRAY <<< "${OUTER_FOLDS}"
declare -A SEEN_HEADS=()
declare -A SEEN_SEEDS=()
declare -A SEEN_FOLDS=()
for head_id in "${HEAD_ARRAY[@]}"; do
  if [[ "${head_id}" != linear && "${head_id}" != two_layer_mlp ]] \
    || [[ -n "${SEEN_HEADS[${head_id}]:-}" ]]; then
    printf 'Head bundle must contain unique admitted common heads.\n' >&2
    exit 1
  fi
  SEEN_HEADS[${head_id}]=1
done
for screen_seed in "${SEED_ARRAY[@]}"; do
  if [[ "${screen_seed}" != 1103 && "${screen_seed}" != 1201 && "${screen_seed}" != 1301 ]] \
    || [[ -n "${SEEN_SEEDS[${screen_seed}]:-}" ]]; then
    printf 'Seed bundle must contain unique frozen screening seeds.\n' >&2
    exit 1
  fi
  SEEN_SEEDS[${screen_seed}]=1
done
for outer_fold in "${FOLD_ARRAY[@]}"; do
  if [[ ! "${outer_fold}" =~ ^[0-4]$ ]] || [[ -n "${SEEN_FOLDS[${outer_fold}]:-}" ]]; then
    printf 'Outer fold bundle must contain unique values in 0..4.\n' >&2
    exit 1
  fi
  SEEN_FOLDS[${outer_fold}]=1
done
for head_id in "${HEAD_ARRAY[@]}"; do
  for screen_seed in "${SEED_ARRAY[@]}"; do
    for outer_fold in "${FOLD_ARRAY[@]}"; do
      "${PYTHON}" scripts/fit_predict_common_cell_heads_study_50000.py shard \
        --embeddings-root "${EMBEDDINGS_ROOT}" \
        --embeddings-relative "${EMBEDDINGS_RELATIVE}" \
        --source "${SOURCE}" \
        --split "${SPLIT}" \
        --output "${STAGE}/shards/${head_id}__seed${screen_seed}__fold${outer_fold}" \
        --expected-embeddings-artifacts-sha256 "${EMBEDDINGS_ARTIFACTS_SHA256}" \
        --expected-source-artifacts-sha256 "${SOURCE_SHA256}" \
        --expected-split-artifacts-sha256 "${SPLIT_SHA256}" \
        --head-id "${head_id}" \
        --screen-seed "${screen_seed}" \
        --outer-fold "${outer_fold}" \
        > "${STAGE}/validation/${head_id}.seed${screen_seed}.fold${outer_fold}.stdout.json" \
        2> "${STAGE}/validation/${head_id}.seed${screen_seed}.fold${outer_fold}.stderr.txt"
    done
  done
done
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/fit_predict_common_cell_heads_study_50000.py \
  tests/unit/test_fit_predict_common_cell_heads_study_50000.py \
  slurm/run_common_cell_head_bundle_50000.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  "${EMBEDDINGS_ROOT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${HEAD_IDS}" "${SCREEN_SEEDS}" "${OUTER_FOLDS}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
heads = sys.argv[2].split()
seeds = [int(value) for value in sys.argv[3].split()]
folds = [int(value) for value in sys.argv[4].split()]
for head in heads:
    for seed in seeds:
        for fold in folds:
            shard = root / f"shards/{head}__seed{seed}__fold{fold}"
            manifest = verify_frozen_tree(shard)
            receipt = json.loads((shard / "prediction_receipt.json").read_text())
            if (
                manifest["metadata"].get("head_id") != head
                or manifest["metadata"].get("screen_seed") != seed
                or manifest["metadata"].get("outer_fold") != fold
                or receipt.get("checkpoint_selection_uses_validation_labels") is not False
                or receipt.get("oof_exact_one_time_coverage") is not True
                or receipt.get("inner_models_used_for_outer_prediction") is not False
                or receipt.get("final_fit_rows") != receipt.get("training_rows")
                or receipt.get("final_fit_donors") != receipt.get("training_donors")
                or receipt.get("final_fit_fixed_epochs") != 30
                or receipt.get("outer_prediction_uses_final_refit") is not True
            ):
                raise SystemExit("common-head shard bundle differs")
freeze_tree(root, {
    "artifact_class": "cell_foundation_common_head_shard_bundle",
    "head_ids": heads,
    "screen_seeds": seeds,
    "outer_folds": folds,
    "logical_shards": len(heads) * len(seeds) * len(folds),
    "metrics_calculated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
