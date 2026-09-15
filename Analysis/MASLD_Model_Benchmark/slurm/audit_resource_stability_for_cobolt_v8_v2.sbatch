#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-audit-210
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-210-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-210-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

REPOSITORY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${REPOSITORY}/Analysis/MASLD_Model_Benchmark
BASELINE=${COBOLT_V8_RESOURCE_BASELINE:?set COBOLT_V8_RESOURCE_BASELINE to the exact frozen baseline}
EXPECTED_PREFIX=${ROOT}/executions/cobolt-v8-resource-baseline-
if [[ "${BASELINE}" != "${EXPECTED_PREFIX}"* ]]; then
  printf 'Invalid Resource baseline binding\n' >&2
  exit 1
fi
TARGET=${ROOT}/executions/cobolt-v8-resource-stability-audit-v2-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Resource stability audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH=src python - "${REPOSITORY}" "${BASELINE}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import time
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.firewall import resource_snapshot

repository = Path(sys.argv[1]).resolve(strict=True)
baseline_path = Path(sys.argv[2]).resolve(strict=True)
stage = Path(sys.argv[3]).resolve(strict=True)
verify_frozen_tree(baseline_path)
baseline = json.loads((baseline_path / "resource_baseline_receipt.json").read_text(encoding="utf-8"))
if (
    baseline.get("schema_version") != "masld-bench-resource-baseline-v1"
    or baseline.get("campaign_id") != "v1-rna-atac-cobolt-5seed-smoke-v8"
    or baseline.get("status") != "pass"
    or baseline.get("outcomes_opened") is not False
):
    raise SystemExit("Resource baseline receipt differs")
started = time.time()
first = resource_snapshot(repository)
latest_mtime = 0.0
latest_path = None
for record in first["files"]:
    path = repository / record["path"]
    stat = path.stat(follow_symlinks=False)
    if stat.st_mtime > latest_mtime:
        latest_mtime = stat.st_mtime
        latest_path = path
second = resource_snapshot(repository)
finished = time.time()
expected_sha256 = baseline["resource_snapshot_sha256"]
if first["snapshot_sha256"] != second["snapshot_sha256"]:
    raise SystemExit("Resource authorities changed during stability double snapshot")
if first["snapshot_sha256"] != expected_sha256:
    raise SystemExit("Resource authorities drifted from the frozen stability baseline")
if latest_path is None:
    raise SystemExit("Resource authority snapshot has no protected files")
if latest_mtime != baseline["latest_protected_file_mtime_unix"]:
    raise SystemExit("latest protected-file timestamp differs from the frozen baseline")
minimum_age_seconds = 1800.0
observed_age_seconds = started - latest_mtime
if observed_age_seconds < minimum_age_seconds:
    raise SystemExit(
        f"Resource stability age is {observed_age_seconds:.3f} seconds; "
        f"requires {minimum_age_seconds:.3f}"
    )
write_json_exclusive(stage / "resource_stability_receipt.json", {
    "schema_version": "masld-bench-resource-stability-audit-v2",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v8",
    "baseline": baseline_path.as_posix(),
    "resource_snapshot_sha256": first["snapshot_sha256"],
    "authority_root_count": len(first["authority_roots"]),
    "authority_file_count": len(first["files"]),
    "latest_protected_file": latest_path.relative_to(repository).as_posix(),
    "latest_protected_file_mtime_unix": latest_mtime,
    "audit_started_unix": started,
    "audit_finished_unix": finished,
    "observed_stability_seconds_at_start": observed_age_seconds,
    "minimum_stability_seconds": minimum_age_seconds,
    "outcomes_opened": False,
    "status": "pass",
})
print(freeze_tree(stage, {
    "artifact_class": "resource_stability_audit",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v8",
    "resource_snapshot_sha256": first["snapshot_sha256"],
    "observed_stability_seconds_at_start": observed_age_seconds,
    "minimum_stability_seconds": minimum_age_seconds,
    "outcomes_opened": False,
    "status": "pass",
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
