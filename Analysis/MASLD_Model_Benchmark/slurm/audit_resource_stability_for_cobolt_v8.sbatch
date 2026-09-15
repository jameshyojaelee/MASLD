#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-audit-208
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-208-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-208-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

REPOSITORY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${REPOSITORY}/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/cobolt-v8-resource-stability-audit-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Resource stability audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH=src python - "${REPOSITORY}" "${STAGE}" <<'PY'
from pathlib import Path
import time
import sys

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.firewall import DEFAULT_RESOURCE_AUTHORITIES, resource_snapshot

repository = Path(sys.argv[1]).resolve(strict=True)
stage = Path(sys.argv[2]).resolve(strict=True)
started = time.time()
snapshot = resource_snapshot(repository)

latest_mtime = 0.0
latest_path = None
directory_count = 0
for relative in DEFAULT_RESOURCE_AUTHORITIES:
    authority = (repository / relative).resolve(strict=True)
    candidates = [authority]
    if authority.is_dir():
        candidates.extend(authority.rglob("*"))
    for path in candidates:
        if path.name == "__pycache__" or path.suffix == ".pyc" or path.is_symlink():
            continue
        stat = path.stat(follow_symlinks=False)
        if path.is_dir():
            directory_count += 1
        if stat.st_mtime > latest_mtime:
            latest_mtime = stat.st_mtime
            latest_path = path
finished = time.time()
minimum_age_seconds = 1800.0
observed_age_seconds = started - latest_mtime
if latest_path is None:
    raise SystemExit("Resource authority tree has no members")
if observed_age_seconds < minimum_age_seconds:
    raise SystemExit(
        f"Resource authority stability gate failed: latest member {latest_path} "
        f"is only {observed_age_seconds:.3f} seconds old"
    )
write_json_exclusive(stage / "resource_stability_receipt.json", {
    "schema_version": "masld-bench-resource-stability-audit-v1",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v8",
    "resource_snapshot_sha256": snapshot["snapshot_sha256"],
    "authority_root_count": len(snapshot["authority_roots"]),
    "authority_file_count": len(snapshot["files"]),
    "authority_directory_count": directory_count,
    "latest_member": latest_path.relative_to(repository).as_posix(),
    "latest_member_mtime_unix": latest_mtime,
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
    "resource_snapshot_sha256": snapshot["snapshot_sha256"],
    "observed_stability_seconds_at_start": observed_age_seconds,
    "minimum_stability_seconds": minimum_age_seconds,
    "outcomes_opened": False,
    "status": "pass",
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
