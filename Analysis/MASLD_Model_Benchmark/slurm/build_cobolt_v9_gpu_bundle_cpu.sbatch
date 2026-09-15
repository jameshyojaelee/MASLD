#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-bundle-212
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-bundle-212-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-bundle-212-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${COBOLT_V9_CANDIDATE:?set COBOLT_V9_CANDIDATE to the exact frozen production candidate}
CANDIDATE_SHA=${COBOLT_V9_CANDIDATE_SHA:?set COBOLT_V9_CANDIDATE_SHA to its ARTIFACTS.json sha256}
STABILITY=${COBOLT_V9_STABILITY:?set COBOLT_V9_STABILITY to the exact stability artifact}
STABILITY_SHA=${COBOLT_V9_STABILITY_SHA:?set COBOLT_V9_STABILITY_SHA to its ARTIFACTS.json sha256}
EXPECTED_PREFIX=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v9--
if [[ "${CANDIDATE}" != "${EXPECTED_PREFIX}"* ]]; then
  printf 'Invalid v9 candidate binding\n' >&2
  exit 1
fi
for digest in "${CANDIDATE_SHA}" "${STABILITY_SHA}"; do
  if [[ "${digest}" == *[!0-9a-f]* || ${#digest} -ne 64 ]]; then
    printf 'Invalid SHA-256 binding\n' >&2
    exit 1
  fi
done

BUILDER=${ROOT}/scripts/build_cobolt_v9_gpu_bundle.py
RUNNER_ROOT=${ROOT}/executions/dispatch-bundles/cobolt-5seed-25run-v2
RUNNER=${RUNNER_ROOT}/run_bundle.py
OUTPUT=${ROOT}/executions/dispatch-bundles/cobolt-5seed-25run-v9-generic5
if [[ -e "${OUTPUT}" ]]; then
  printf 'Refusing existing bundle output: %s\n' "${OUTPUT}" >&2
  exit 1
fi

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/cobolt-v9-bundle-${SLURM_JOB_ID}
cd "${ROOT}"
printf '%s  %s\n' "${CANDIDATE_SHA}" "${CANDIDATE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${STABILITY_SHA}" "${STABILITY}/ARTIFACTS.json" | sha256sum --check --strict
PYTHONPATH=src python - "${STABILITY}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(root)
receipt = json.loads((root / "resource_stability_receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("schema_version") != "masld-bench-cobolt-v9-resource-stability-v1"
    or receipt.get("campaign_id") != "v1-rna-atac-cobolt-5seed-smoke-v9"
    or receipt.get("status") != "pass"
    or receipt.get("outcomes_opened") is not False
):
    raise SystemExit("Cobolt v9 stability receipt differs")
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
  --runner "${RUNNER}"
PYTHONPATH=src python "${RUNNER}" --manifest "${OUTPUT}/bundle_manifest.json" --preflight
sha256sum "${OUTPUT}/ARTIFACTS.json" "${OUTPUT}/model_work_212.sbatch"
