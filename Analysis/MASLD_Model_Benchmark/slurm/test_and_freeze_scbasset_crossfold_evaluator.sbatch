#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-039
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-039-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-039-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
EVALUATOR=${ROOT}/scripts/evaluate_scbasset_valid_crossfold.py
LOCK_BUILDER=${ROOT}/scripts/freeze_scbasset_crossfold_prediction_lock.py
BASE_EVALUATOR=${ROOT}/scripts/evaluate_scbasset_valid.py
TEST=${ROOT}/tests/unit/test_evaluate_scbasset_valid_crossfold.py
LOCK_TEST=${ROOT}/tests/unit/test_freeze_scbasset_crossfold_prediction_lock.py
BUNDLE=${ROOT}/executions/dispatch-bundles/scbasset-fold1to4-chain-v2-generic1
PARENT=${ROOT}/executions/dispatch-bundles/scbasset-fold1to4-chain-v1
BASE=${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230
IMAGE=${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif
CORE=${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507
TARGET=${ROOT}/executions/scbasset-crossfold-evaluator-source-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite scBasset evaluator source closure: %s\n' "${TARGET}" >&2
  exit 1
fi

declare -a INPUTS=(
  "${ROOT}/executions/scbasset-donor1_genomic1-inputs-21067171"
  "${ROOT}/executions/scbasset-donor2_genomic2-inputs-21067172"
  "${ROOT}/executions/scbasset-donor3_genomic3-inputs-21067173"
  "${ROOT}/executions/scbasset-donor4_genomic4-inputs-21067174"
)
declare -a INPUT_HASHES=(
  5a608a40763e0d72e0af1a284690c762f54379a8230e06889361176935bb88c3
  1b449406fba08190a94bb2750eaf558e7014b187fb239c998659df68561f610f
  bb9b1ce8f4e84ddd4509e97ee593ddca5a39e4b66556aa228c057c50a35f9678
  4796aa202ec5f2a171f7422fcb80544d66f9e288372816a3e856cf688e8d4505
)

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH=${ROOT}:${ROOT}/src python -m py_compile \
  "${EVALUATOR}" "${LOCK_BUILDER}" "${TEST}" "${LOCK_TEST}"
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest \
  tests.unit.test_freeze_scbasset_crossfold_prediction_lock
printf '%s  %s\n' \
  b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 "${BASE}/ARTIFACTS.json" \
  96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa "${IMAGE}" \
  61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 "${CORE}/ARTIFACTS.json" \
  | sha256sum --check --strict
module purge
module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:ro" \
  --env "PYTHONPATH=${ROOT}:${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m unittest \
  tests.unit.test_evaluate_scbasset_valid_crossfold \
  tests.unit.test_evaluate_scbasset_valid
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}:${ROOT}/src python - \
  "${BUNDLE}" "${PARENT}" "${INPUTS[@]}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

bundle = Path(sys.argv[1]).resolve(strict=True)
parent = Path(sys.argv[2]).resolve(strict=True)
inputs = [Path(value).resolve(strict=True) for value in sys.argv[3:]]
verify_frozen_tree(bundle)
verify_frozen_tree(parent)
for fold, root in enumerate(inputs, start=1):
    verify_frozen_tree(root)
    summary = json.loads((root / "inputs/summary.json").read_text(encoding="utf-8"))
    if (
        summary["status"] != "pass"
        or summary["donor_test_fold"] != fold
        or summary["donor_valid_fold"] != (fold + 1) % 5
        or summary["held_donor_atac_exported"]
        or summary["genomic_test_atac_exported"]
    ):
        raise SystemExit(f"fold {fold} evaluator input contract differs")
PY
for index in 0 1 2 3; do
  printf '%s  %s\n' "${INPUT_HASHES[${index}]}" "${INPUTS[${index}]}/ARTIFACTS.json" \
    | sha256sum --check --strict
done

install -d -m 0750 "${STAGE}/source"
install -m 0640 "${EVALUATOR}" "${STAGE}/source/evaluate_scbasset_valid_crossfold.py"
install -m 0640 "${LOCK_BUILDER}" "${STAGE}/source/freeze_scbasset_crossfold_prediction_lock.py"
install -m 0640 "${BASE_EVALUATOR}" "${STAGE}/source/evaluate_scbasset_valid.py"
install -m 0640 "${TEST}" "${STAGE}/source/test_evaluate_scbasset_valid_crossfold.py"
install -m 0640 "${LOCK_TEST}" "${STAGE}/source/test_freeze_scbasset_crossfold_prediction_lock.py"
PYTHONPATH=${ROOT}/src python - \
  "${STAGE}" "${EVALUATOR}" "${LOCK_BUILDER}" "${BASE_EVALUATOR}" "${TEST}" "${LOCK_TEST}" \
  "${BUNDLE}" "${PARENT}" "${INPUTS[@]}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.hashing import sha256_file

stage = Path(sys.argv[1]).resolve(strict=True)
evaluator = Path(sys.argv[2]).resolve(strict=True)
lock_builder = Path(sys.argv[3]).resolve(strict=True)
base_evaluator = Path(sys.argv[4]).resolve(strict=True)
test = Path(sys.argv[5]).resolve(strict=True)
lock_test = Path(sys.argv[6]).resolve(strict=True)
bundle = Path(sys.argv[7]).resolve(strict=True)
parent = Path(sys.argv[8]).resolve(strict=True)
inputs = [Path(value).resolve(strict=True) for value in sys.argv[9:]]
records = {
    "evaluator": {"path": evaluator.as_posix(), "sha256": sha256_file(evaluator)},
    "prediction_lock_builder": {"path": lock_builder.as_posix(), "sha256": sha256_file(lock_builder)},
    "base_evaluator": {"path": base_evaluator.as_posix(), "sha256": sha256_file(base_evaluator)},
    "unit_test": {"path": test.as_posix(), "sha256": sha256_file(test)},
    "lock_unit_test": {"path": lock_test.as_posix(), "sha256": sha256_file(lock_test)},
    "active_bundle": {"path": bundle.as_posix(), "sha256": sha256_file(bundle / "ARTIFACTS.json")},
    "parent_bundle": {"path": parent.as_posix(), "sha256": sha256_file(parent / "ARTIFACTS.json")},
    "inputs": [
        {"path": path.as_posix(), "sha256": sha256_file(path / "ARTIFACTS.json")}
        for path in inputs
    ],
}
write_json_exclusive(stage / "source_lock.json", {
    "schema_version": "masld-bench-scbasset-crossfold-evaluator-source-lock-v1",
    "model_id": "scbasset",
    "dataset_id": "gse296875",
    "seed": 20260824,
    "folds": [1, 2, 3, 4],
    "records": records,
    "predictions_available_at_freeze": False,
    "outcomes_read": False,
    "champion_claim_allowed": False,
})
print(freeze_tree(stage, {
    "artifact_class": "scbasset_crossfold_evaluator_source_lock",
    "model_id": "scbasset",
    "dataset_id": "gse296875",
    "seed": 20260824,
    "folds": [1, 2, 3, 4],
    "evaluator_sha256": records["evaluator"]["sha256"],
    "prediction_lock_builder_sha256": records["prediction_lock_builder"]["sha256"],
    "outcomes_read": False,
    "status": "pass",
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${EVALUATOR}" "${LOCK_BUILDER}" "${TEST}"
