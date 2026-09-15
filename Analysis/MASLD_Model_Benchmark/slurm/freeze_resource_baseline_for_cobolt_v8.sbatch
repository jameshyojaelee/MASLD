#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-audit-209
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-209-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-209-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

REPOSITORY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${REPOSITORY}/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/cobolt-v8-resource-baseline-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Resource baseline: %s\n' "${TARGET}" >&2
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
from masld_bench.firewall import resource_snapshot

repository = Path(sys.argv[1]).resolve(strict=True)
stage = Path(sys.argv[2]).resolve(strict=True)
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
if first["snapshot_sha256"] != second["snapshot_sha256"]:
    raise SystemExit("Resource authorities changed during baseline double snapshot")
if latest_path is None:
    raise SystemExit("Resource authority snapshot has no protected files")
write_json_exclusive(stage / "resource_baseline_receipt.json", {
    "schema_version": "masld-bench-resource-baseline-v1",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v8",
    "resource_snapshot_sha256": first["snapshot_sha256"],
    "authority_root_count": len(first["authority_roots"]),
    "authority_file_count": len(first["files"]),
    "latest_protected_file": latest_path.relative_to(repository).as_posix(),
    "latest_protected_file_mtime_unix": latest_mtime,
    "baseline_started_unix": started,
    "baseline_finished_unix": finished,
    "outcomes_opened": False,
    "status": "pass",
})
print(freeze_tree(stage, {
    "artifact_class": "resource_stability_baseline",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v8",
    "resource_snapshot_sha256": first["snapshot_sha256"],
    "latest_protected_file_mtime_unix": latest_mtime,
    "outcomes_opened": False,
    "status": "pass",
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
