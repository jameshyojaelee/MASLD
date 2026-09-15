#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-bundle-358
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-bundle-358-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-bundle-358-%j.err

# Build and freeze the v12 B6K GPU dispatch bundle for the 25 Cobolt fits.
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
V12=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v12--3972f43112b8ecc4
V12_SHA=9b878eb08b3fc8a6012f473694b778a771da1bdd976f59e2966776f0bda3519f
V11=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v11--e17be3dab9ba48cb
V11_SHA=dc4f962c841f6daaf50eadcbd476d1ed26eb0a57deada67222fa3ff998c71f47
SCOPE=${ROOT}/executions/cobolt-v12-revision-scope-validation-21111670
RUNNER=${ROOT}/executions/dispatch-bundles/cobolt-5seed-25run-v2/run_bundle.py
BUILDER=${ROOT}/scripts/build_cobolt_v12_gpu_bundle.py
OUTPUT=${ROOT}/executions/dispatch-bundles/cobolt-5seed-25run-v12-b6k

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
printf '%s  %s\n' "${V12_SHA}" "${V12}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${V11_SHA}" "${V11}/ARTIFACTS.json" | sha256sum --check --strict

# The v12 scope validation must have passed AND its separation must still be current.
PYTHONPATH=src python - "${SCOPE}" <<'PY'
import json, sys
from pathlib import Path
from masld_bench.artifacts import verify_frozen_tree
scope = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(scope)
r = json.loads((scope / "receipt.json").read_text())
if r.get("status") != "passed_firewall_baseline_only_revision":
    raise SystemExit(f"v12 scope validation did not pass: {r.get('status')}")
if r.get("resource_firewall_still_current") is not True:
    raise SystemExit("v12 firewall baseline was already stale at scope validation")
if r.get("scientific_content_identical_to_v11") is not True:
    raise SystemExit("v12 scientific content is not identical to v11")
print("v12 scope validation verified")
PY

PYTHONPATH=src python "${BUILDER}" \
  --candidate "${V12}" \
  --superseded-candidate "${V11}" \
  --output "${OUTPUT}" \
  --runner "${RUNNER}"
