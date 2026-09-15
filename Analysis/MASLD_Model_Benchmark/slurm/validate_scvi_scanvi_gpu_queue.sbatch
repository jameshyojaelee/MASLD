#!/usr/bin/env bash
#SBATCH --job-name=model-queue-148
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/scvi-tools-1.4.2-inductive-portable-21079606
RUNTIME_SHA256=d885ada24fdad9f6e3d57a92ff870ca8ff05f07701ac6b0762f55a8cb19db61a
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
PREFLIGHT=${ROOT}/executions/model-preflight-146-21079904
PREFLIGHT_SHA256=cb4eae65816d03b2877afb0fbd8f802f29c252feacc5f98916c79ef64ce32267
STAGED_QUEUE=${ROOT}/config/campaigns/gpu_bundle_queue_staging/055-scvi-scanvi-50k-study-held-3seed-v1.json
STAGED_QUEUE_SHA256=8035c8720f145748cd885887e23d35d74768f4deb340f25c33c441ab16752a0c
LIVE_QUEUE=${ROOT}/config/campaigns/gpu_bundle_queue/055-scvi-scanvi-50k-study-held-3seed-v1.json
WRAPPER=${ROOT}/slurm/run_scvi_scanvi_study_50000_seed_bundle.sbatch
WRAPPER_SHA256=c4661c4d4313419bad9e96f1b06bdbeae3779e7600c0b18c8422f2b4620f5dbe
CONTRACT=${ROOT}/config/artifacts/models/cell_bulk_classical_baselines/scvi_scanvi_1_4_2_frozen_screen_contract.json
CONTRACT_SHA256=32d86487d5ffc61a8d2f74d49da7eb9ba09fff67e9448506e52c8783f6bc0c04
ACTIVATION=${ROOT}/executions/scvi-scanvi-gpu-queue-validation-v1
CLAIM=${ROOT}/executions/gpu-bundle-dispatch-state/claims/scvi-scanvi-50k-study-held-3seed-v1.json
TARGET=${ACTIVATION}
STAGE=${TARGET}.staging-${SLURM_JOB_ID}
FAILED=${TARGET}.failed-${SLURM_JOB_ID}

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite scVI/scANVI GPU queue validation.\n' >&2
  exit 1
fi
if [[ -e "${LIVE_QUEUE}" || -e "${CLAIM}" ]]; then
  printf 'scVI/scANVI queue or dispatcher claim exists before validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"

printf '%s  %s\n' \
  "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  "${PREFLIGHT_SHA256}" "${PREFLIGHT}/ARTIFACTS.json" \
  "${STAGED_QUEUE_SHA256}" "${STAGED_QUEUE}" \
  "${WRAPPER_SHA256}" "${WRAPPER}" \
  "${CONTRACT_SHA256}" "${CONTRACT}" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"
bash -n "${WRAPPER}"

module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH=${RUNTIME}/site-packages:${ROOT}:${ROOT}/src
python -m unittest -v \
  tests.unit.test_gpu_bundle_dispatcher \
  tests.unit.test_scvi_scanvi_inductive \
  tests.unit.test_fit_predict_scvi_scanvi_study_50000 \
  tests.unit.test_aggregate_scvi_scanvi_study_50000 \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

python - "${STAGED_QUEUE}" "${ROOT}" "${LIVE_QUEUE}" "${CLAIM}" \
  "${ACTIVATION}/ARTIFACTS.json" \
  > "${STAGE}/validation/prepublication.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree
from scripts.gpu_bundle_dispatcher import load_item

queue = Path(sys.argv[1])
root = Path(sys.argv[2]).resolve(strict=True)
live_queue = Path(sys.argv[3])
claim = Path(sys.argv[4])
activation_artifacts = Path(sys.argv[5])
item = load_item(queue, root)
wrapper = item["_wrapper"].read_text(encoding="utf-8")
fit_script = (root / "scripts/fit_predict_scvi_scanvi_study_50000.py").read_text(
    encoding="utf-8"
)
required_header = {
    "#SBATCH --job-name=model-work-151",
    "#SBATCH --partition=gpu",
    "#SBATCH --account=nslab",
    "#SBATCH --qos=nslab",
    "#SBATCH --gres=gpu:l40s:1",
    "#SBATCH --cpus-per-task=8",
    "#SBATCH --mem=96G",
    "#SBATCH --time=72:00:00",
}
expected_paths = {
    "executions/environments/scvi-tools-1.4.2-inductive-portable-21079606/ARTIFACTS.json",
    "executions/model-data-050-20260824T211318Z-21077049/fixture/ARTIFACTS.json",
    "executions/model-data-052-21077595/split/ARTIFACTS.json",
    "executions/model-preflight-146-21079904/ARTIFACTS.json",
    "config/artifacts/models/cell_bulk_classical_baselines/scvi_scanvi_1_4_2_frozen_screen_contract.json",
    "executions/scvi-scanvi-gpu-queue-validation-v1/ARTIFACTS.json",
}
existing_required = [
    (root / relative).resolve(strict=False)
    for relative in item["required_paths"]
    if relative != "executions/scvi-scanvi-gpu-queue-validation-v1/ARTIFACTS.json"
]
if (
    item["_ready"]
    or item["bundle_id"] != "scvi-scanvi-50k-study-held-3seed-v1"
    or item["priority"] != 55
    or item["logical_tasks"] != 15
    or item["family"] != "cell_transcriptome"
    or not item["enabled"]
    or set(item["required_paths"]) != expected_paths
    or not all(path.is_file() for path in existing_required)
    or activation_artifacts.exists()
    or live_queue.exists()
    or claim.exists()
    or not required_header <= set(wrapper.splitlines())
    or "#SBATCH --array" in wrapper
    or "--qos=innovation" in wrapper
    or any(line.strip().startswith("sbatch ") for line in wrapper.splitlines())
    or "for SEED in 20260824 20260825 20260826; do" not in wrapper
    or wrapper.count("--seed-root \"${STAGE}/seed-") != 3
    or "aggregate_scvi_scanvi_study_50000.py" not in wrapper
):
    raise SystemExit("staged scVI/scANVI queue contract differs")
for forbidden in (
    "load_query_data",
    "prepare_query_anndata",
    "scarches",
    "SCVI.load_query_data",
    "SCANVI.load_query_data",
):
    if forbidden in fit_script:
        raise SystemExit(f"forbidden query adaptation API present: {forbidden}")
for frozen_root in (
    root / "executions/environments/scvi-tools-1.4.2-inductive-portable-21079606",
    root / "executions/model-data-050-20260824T211318Z-21077049/fixture",
    root / "executions/model-data-052-21077595/split",
    root / "executions/model-preflight-146-21079904",
):
    verify_frozen_tree(frozen_root)
print(
    json.dumps(
        {
            "bundle_id": item["bundle_id"],
            "ready_before_validation": False,
            "only_activation_gate_missing": True,
            "live_queue_absent": True,
            "dispatcher_claim_absent": True,
            "gpu_submitted": False,
            "project_data_fit": False,
            "outcomes_read": False,
        },
        sort_keys=True,
    )
)
PY
cp "${RUNTIME}/validation/pip-freeze.txt" \
  "${STAGE}/validation/runtime-pip-freeze.txt"
sha256sum \
  scripts/gpu_bundle_dispatcher.py \
  scripts/fit_predict_scvi_scanvi_study_50000.py \
  scripts/aggregate_scvi_scanvi_study_50000.py \
  src/masld_bench/adapters/scvi_scanvi_inductive.py \
  tests/unit/test_gpu_bundle_dispatcher.py \
  tests/unit/test_scvi_scanvi_inductive.py \
  tests/unit/test_fit_predict_scvi_scanvi_study_50000.py \
  tests/unit/test_aggregate_scvi_scanvi_study_50000.py \
  slurm/run_scvi_scanvi_study_50000_seed_bundle.sbatch \
  slurm/validate_scvi_scanvi_gpu_queue.sbatch \
  config/campaigns/gpu_bundle_queue_staging/055-scvi-scanvi-50k-study-held-3seed-v1.json \
  config/artifacts/models/cell_bulk_classical_baselines/scvi_scanvi_1_4_2_frozen_screen_contract.json \
  > "${STAGE}/validation/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${STAGE}" "${STAGED_QUEUE_SHA256}" \
  "${WRAPPER_SHA256}" "${RUNTIME_SHA256}" "${PREFLIGHT_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
receipt = {
    "schema_version": "masld-bench-scvi-scanvi-gpu-queue-validation-v1",
    "status": "passed",
    "bundle_id": "scvi-scanvi-50k-study-held-3seed-v1",
    "priority": 55,
    "logical_fold_seed_fits": 15,
    "seeds": [20260824, 20260825, 20260826],
    "models_per_fit": ["scvi_baseline", "scanvi_baseline"],
    "staged_queue_sha256": sys.argv[2],
    "wrapper_sha256": sys.argv[3],
    "runtime_artifacts_sha256": sys.argv[4],
    "preflight_artifacts_sha256": sys.argv[5],
    "partition": "gpu",
    "account": "nslab",
    "qos": "nslab",
    "gres": "gpu:l40s:1",
    "cpus": 8,
    "memory": "96G",
    "time": "72:00:00",
    "array": False,
    "job_name": "model-work-151",
    "dispatcher_only": True,
    "shared_gate": {"cap": 5, "max_running": 4, "max_pending": 1},
    "direct_held_study_query_only": True,
    "query_training_allowed": False,
    "query_adaptation_allowed": False,
    "common_head_fit": False,
    "metrics_run": False,
    "evaluation_labels_read": False,
    "sealed_outcomes_read": False,
    "gpu_submitted": False,
    "live_queue_mutated": False,
}
(stage / "validation/receipt.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
freeze_tree(
    stage,
    {
        "artifact_class": "scvi_scanvi_gpu_queue_contract_validation",
        "bundle_id": receipt["bundle_id"],
        "priority": receipt["priority"],
        "staged_queue_sha256": sys.argv[2],
        "wrapper_sha256": sys.argv[3],
        "gpu_submitted": False,
        "live_queue_mutated": False,
        "status": "passed",
    },
)
PY
PYTHONPATH=${ROOT}/src python - "${STAGE}" "${TARGET}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import publish_directory_noreplace

publish_directory_noreplace(Path(sys.argv[1]), Path(sys.argv[2]))
PY
PYTHONPATH=${ROOT}:${ROOT}/src python - "${STAGED_QUEUE}" "${ROOT}" \
  "${TARGET}" "${LIVE_QUEUE}" "${CLAIM}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree
from scripts.gpu_bundle_dispatcher import load_item

verify_frozen_tree(Path(sys.argv[3]))
item = load_item(Path(sys.argv[1]), Path(sys.argv[2]).resolve(strict=True))
if (
    not item["_ready"]
    or item["bundle_id"] != "scvi-scanvi-50k-study-held-3seed-v1"
    or Path(sys.argv[4]).exists()
    or Path(sys.argv[5]).exists()
):
    raise SystemExit("scVI/scANVI staged queue is not ready after validation")
print(item["_queue_sha256"])
PY

trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json"
