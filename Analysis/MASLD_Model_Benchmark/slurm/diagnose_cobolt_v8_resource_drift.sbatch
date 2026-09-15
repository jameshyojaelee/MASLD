#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-audit-211
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-211-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-211-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

REPOSITORY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${REPOSITORY}/Analysis/MASLD_Model_Benchmark
BASELINE=${ROOT}/executions/cobolt-v8-resource-baseline-21075101
TARGET=${ROOT}/executions/cobolt-v8-resource-drift-diagnostic-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Resource drift diagnostic: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH=src python - "${REPOSITORY}" "${BASELINE}" "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.firewall import DEFAULT_RESOURCE_AUTHORITIES

repository = Path(sys.argv[1]).resolve(strict=True)
baseline_path = Path(sys.argv[2]).resolve(strict=True)
stage = Path(sys.argv[3]).resolve(strict=True)
verify_frozen_tree(baseline_path)
baseline = json.loads((baseline_path / "resource_baseline_receipt.json").read_text(encoding="utf-8"))
baseline_started = float(baseline["baseline_started_unix"])
baseline_finished = float(baseline["baseline_finished_unix"])
audit_started = 1787603641.0  # 2026-08-24T16:34:01-04:00 from SLURM accounting.

records = []
for relative in DEFAULT_RESOURCE_AUTHORITIES:
    authority = (repository / relative).resolve(strict=True)
    pending = [authority]
    while pending:
        path = pending.pop()
        if path.name == "__pycache__" or path.suffix == ".pyc" or path.is_symlink():
            continue
        stat = path.stat(follow_symlinks=False)
        kind = "directory" if path.is_dir() else "file"
        if stat.st_mtime >= baseline_started:
            records.append({
                "path": path.relative_to(repository).as_posix(),
                "kind": kind,
                "mtime_unix": stat.st_mtime,
                "modified_after_baseline_started": stat.st_mtime >= baseline_started,
                "modified_after_baseline_finished": stat.st_mtime >= baseline_finished,
                "modified_after_audit_started": stat.st_mtime >= audit_started,
            })
        if path.is_dir():
            pending.extend(sorted(path.iterdir(), key=lambda item: item.name, reverse=True))
records.sort(key=lambda item: (-float(item["mtime_unix"]), str(item["path"])))
fields = (
    "path", "kind", "mtime_unix", "modified_after_baseline_started",
    "modified_after_baseline_finished", "modified_after_audit_started",
)
with (stage / "protected_members_modified_since_baseline.tsv").open(
    "x", encoding="utf-8", newline=""
) as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
files_after_audit = [
    record for record in records
    if record["kind"] == "file" and record["modified_after_audit_started"]
]
directories_after_audit = [
    record for record in records
    if record["kind"] == "directory" and record["modified_after_audit_started"]
]
write_json_exclusive(stage / "diagnostic_receipt.json", {
    "schema_version": "masld-bench-resource-drift-diagnostic-v1",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v8",
    "baseline": baseline_path.as_posix(),
    "baseline_snapshot_sha256": baseline["resource_snapshot_sha256"],
    "baseline_started_unix": baseline_started,
    "baseline_finished_unix": baseline_finished,
    "audit_started_unix": audit_started,
    "protected_members_modified_after_baseline_started": len(records),
    "protected_files_modified_after_audit_started": len(files_after_audit),
    "protected_directories_modified_after_audit_started": len(directories_after_audit),
    "newest_members": records[:20],
    "outcomes_opened": False,
    "authorities_modified": False,
    "status": "diagnostic_complete",
})
print(freeze_tree(stage, {
    "artifact_class": "resource_drift_diagnostic",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v8",
    "protected_files_modified_after_audit_started": len(files_after_audit),
    "protected_directories_modified_after_audit_started": len(directories_after_audit),
    "outcomes_opened": False,
    "authorities_modified": False,
    "status": "diagnostic_complete",
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
