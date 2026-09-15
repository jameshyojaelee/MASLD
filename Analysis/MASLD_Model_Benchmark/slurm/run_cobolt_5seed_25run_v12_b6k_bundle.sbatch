#!/usr/bin/env bash
#SBATCH --job-name=model-work-303
#SBATCH --account=nslab
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --gres=gpu:b6k:1
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-work-303-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-work-303-%j.err

# The 25 Cobolt GPU fits for campaign v12, on B6K.
#
# One job, one B6K device, all 25 fits serially through the same reviewed
# run_bundle.py driver the L40S bundles used. Hardware homogeneity within the
# rectangle therefore holds by construction.
#
# Direct sbatch: the GPU dependency dispatcher hardcodes gpu:l40s:1 in four
# places and cannot route B6K. Job ids from this wrapper never enter its roster,
# so job 21100574 and the four inherited lanes are unaffected.
#
# Runtime gpu_rna_atac_torch_b6k, validated on sm_120 in job 21107495: forward/
# backward equation, query invariance and tensor state restore re-established
# from scratch, plus cross-process determinism. The L40S validation is NOT
# reused. torch is the only delta from the validated L40S runtime.

set -euo pipefail
umask 027
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BUNDLE=${ROOT}/executions/dispatch-bundles/cobolt-5seed-25run-v12-b6k
CANDIDATE=${COBOLT_V12_CANDIDATE:?set COBOLT_V12_CANDIDATE to the frozen v12 candidate}
RUNNER=${ROOT}/executions/dispatch-bundles/cobolt-5seed-25run-v2/run_bundle.py
RUNNER_SHA=b710857a049e25e88a2d653566146c6c352ebef67c983e8fefac52cae7f64274

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH=src python - "${ROOT}" "${CANDIDATE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.planner import load_frozen_plan, source_lock
root = Path(sys.argv[1]).resolve(strict=True)
plan = load_frozen_plan(Path(sys.argv[2]))
if plan["campaign"]["campaign_id"] != "v1-rna-atac-cobolt-5seed-smoke-v12":
    raise SystemExit("candidate is not the v12 campaign; zero runs started")
included_paths = plan["source_lock"].get("included_paths")
if source_lock(root, included_paths=included_paths) != plan["source_lock"]:
    raise SystemExit("source tree differs from frozen v12 candidate; zero runs started")
PY
printf '%s  %s\n' "${RUNNER_SHA}" "${RUNNER}" | sha256sum --check --strict
PYTHONPATH=src python "${RUNNER}" \
  --manifest "${BUNDLE}/bundle_manifest.json" \
  --execute --execution-id "${SLURM_JOB_ID}"
