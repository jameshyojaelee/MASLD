#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-audit-211
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-211-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-211-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

REPOSITORY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${REPOSITORY}/Analysis/MASLD_Model_Benchmark
PREFLIGHT=${COBOLT_V9_PREFLIGHT:?set COBOLT_V9_PREFLIGHT to the exact frozen v9 preflight}
EXPECTED_PREFIX=${ROOT}/executions/cobolt-v9-campaign-preflight-
if [[ "${PREFLIGHT}" != "${EXPECTED_PREFIX}"* ]]; then
  printf 'Invalid v9 preflight binding\n' >&2
  exit 1
fi
TARGET=${ROOT}/executions/cobolt-v9-resource-stability-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite v9 Resource stability audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH=src python - "${REPOSITORY}" "${PREFLIGHT}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
import time

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.firewall import resource_snapshot
from masld_bench.planner import load_frozen_plan

repository = Path(sys.argv[1]).resolve(strict=True)
preflight = Path(sys.argv[2]).resolve(strict=True)
stage = Path(sys.argv[3]).resolve(strict=True)
verify_frozen_tree(preflight)
receipt = json.loads((preflight / "preflight_receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("schema_version") != "masld-bench-cobolt-v9-campaign-preflight-v1"
    or receipt.get("campaign_id") != "v1-rna-atac-cobolt-5seed-smoke-v9"
    or receipt.get("status") != "pass"
    or receipt.get("outcomes_opened") is not False
):
    raise SystemExit("v9 preflight receipt differs")
candidate = preflight / receipt["preflight_candidate"]
verify_frozen_tree(candidate)
plan = load_frozen_plan(candidate)
expected = plan["resource_firewall"]["snapshot_sha256"]
if expected != receipt["resource_snapshot_sha256"]:
    raise SystemExit("preflight Resource snapshot bindings disagree")
started = time.time()
first = resource_snapshot(repository)
second = resource_snapshot(repository)
if first["snapshot_sha256"] != second["snapshot_sha256"]:
    raise SystemExit("Resource authorities changed during stability double snapshot")
if first["snapshot_sha256"] != expected:
    raise SystemExit("Resource authorities drifted after the v9 preflight")
latest_mtime = 0.0
latest_path = None
for record in first["files"]:
    path = repository / record["path"]
    observed = path.stat(follow_symlinks=False).st_mtime
    if observed > latest_mtime:
        latest_mtime = observed
        latest_path = path
if latest_path is None:
    raise SystemExit("Resource authority snapshot has no protected files")
minimum_age_seconds = 1800.0
observed_age_seconds = started - latest_mtime
if observed_age_seconds < minimum_age_seconds:
    raise SystemExit(
        f"Resource stability age is {observed_age_seconds:.3f} seconds; "
        f"requires {minimum_age_seconds:.3f}"
    )
write_json_exclusive(stage / "resource_stability_receipt.json", {
    "schema_version": "masld-bench-cobolt-v9-resource-stability-v1",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v9",
    "preflight": preflight.as_posix(),
    "preflight_artifacts_sha256": receipt["candidate_artifacts_sha256"],
    "resource_snapshot_sha256": expected,
    "authority_root_count": len(first["authority_roots"]),
    "authority_file_count": len(first["files"]),
    "latest_protected_file": latest_path.relative_to(repository).as_posix(),
    "latest_protected_file_mtime_unix": latest_mtime,
    "observed_stability_seconds_at_start": observed_age_seconds,
    "minimum_stability_seconds": minimum_age_seconds,
    "outcomes_opened": False,
    "status": "pass",
})
print(freeze_tree(stage, {
    "artifact_class": "cobolt_v9_resource_stability",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v9",
    "resource_snapshot_sha256": expected,
    "observed_stability_seconds_at_start": observed_age_seconds,
    "minimum_stability_seconds": minimum_age_seconds,
    "outcomes_opened": False,
    "status": "pass",
}))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
