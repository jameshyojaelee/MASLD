#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-bundle-215
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-bundle-215-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-bundle-215-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v9--edf243d20e2b8f58
CANDIDATE_SHA=ac24e1a8da4fde1c35aafa1370659155eded36811099b40f3ee39a4fc2866862
STABILITY=${ROOT}/executions/cobolt-v9-frozen-candidate-stability-21087850
STABILITY_SHA=b8b6003576becb93b530c20465d9022ee334bae5bbdcb660ba340836d0831233
BUILDER=${ROOT}/scripts/build_cobolt_v9_gpu_bundle.py
BUILDER_SHA=837405a39762ba824770f08c45c610b419f5e6b490ceec69f7b21f6836998cee
RUNNER_ROOT=${ROOT}/executions/dispatch-bundles/cobolt-5seed-25run-v2
RUNNER=${RUNNER_ROOT}/run_bundle.py
RUNNER_SHA=b710857a049e25e88a2d653566146c6c352ebef67c983e8fefac52cae7f64274
RUNNER_ROOT_SHA=4f9695fc7a7397ca5437ba6235e0d5f178599c3cbbe0e41b4a7190a949ec8a7a
BUNDLE_ID=cobolt-5seed-25run-v9-generic6
OUTPUT=${ROOT}/executions/dispatch-bundles/${BUNDLE_ID}

if [[ -e ${OUTPUT} ]]; then
  printf 'Refusing existing bundle output: %s\n' "${OUTPUT}" >&2
  exit 1
fi
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/cobolt-v9-generic6-${SLURM_JOB_ID}
cd "${ROOT}"
printf '%s  %s\n' "${CANDIDATE_SHA}" "${CANDIDATE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${STABILITY_SHA}" "${STABILITY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${BUILDER_SHA}" "${BUILDER}" | sha256sum --check --strict
printf '%s  %s\n' "${RUNNER_SHA}" "${RUNNER}" | sha256sum --check --strict
printf '%s  %s\n' "${RUNNER_ROOT_SHA}" "${RUNNER_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
PYTHONPATH=src python - "${CANDIDATE}" "${STABILITY}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree
from masld_bench.hashing import sha256_file

candidate, stability = map(lambda value: Path(value).resolve(strict=True), sys.argv[1:])
verify_frozen_tree(candidate)
verify_frozen_tree(stability)
receipt = json.loads((stability / "resource_stability_receipt.json").read_text())
if (
    receipt.get("status") != "pass"
    or receipt.get("candidate_artifacts_sha256") != sha256_file(candidate / "ARTIFACTS.json")
    or receipt.get("outcomes_opened") is not False
):
    raise SystemExit("Cobolt v9 stability binding differs")
PY
PYTHONPATH=src python -m py_compile "${RUNNER}" "${BUILDER}"
PYTHONPATH=src python "${RUNNER}" --self-test
PYTHONPATH=src python - "${RUNNER_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
PY
PYTHONPATH=src python "${BUILDER}" \
  --candidate "${CANDIDATE}" \
  --output "${OUTPUT}" \
  --runner "${RUNNER}" \
  --bundle-id "${BUNDLE_ID}"
PYTHONPATH=src python "${RUNNER}" --manifest "${OUTPUT}/bundle_manifest.json" --preflight
sha256sum "${OUTPUT}/ARTIFACTS.json" "${OUTPUT}/model_work_212.sbatch"
