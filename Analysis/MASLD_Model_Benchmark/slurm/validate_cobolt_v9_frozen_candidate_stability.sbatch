#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-audit-212
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-212-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-audit-212-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

REPOSITORY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${REPOSITORY}/Analysis/MASLD_Model_Benchmark
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v9--edf243d20e2b8f58
CANDIDATE_SHA=ac24e1a8da4fde1c35aafa1370659155eded36811099b40f3ee39a4fc2866862
TARGET=${ROOT}/executions/cobolt-v9-frozen-candidate-stability-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite frozen-candidate stability audit.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
printf '%s  %s\n' "${CANDIDATE_SHA}" "${CANDIDATE}/ARTIFACTS.json" | sha256sum --check --strict
PYTHONPATH=src python - "${REPOSITORY}" "${CANDIDATE}" "${CANDIDATE_SHA}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
import time

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.firewall import resource_snapshot
from masld_bench.hashing import sha256_file
from masld_bench.planner import load_frozen_plan, source_lock

repository = Path(sys.argv[1]).resolve(strict=True)
candidate = Path(sys.argv[2]).resolve(strict=True)
candidate_sha = sys.argv[3]
stage = Path(sys.argv[4]).resolve(strict=True)
verify_frozen_tree(candidate)
if sha256_file(candidate / "ARTIFACTS.json") != candidate_sha:
    raise SystemExit("frozen candidate identity differs")
plan = load_frozen_plan(candidate)
if (
    plan["campaign"]["campaign_id"] != "v1-rna-atac-cobolt-5seed-smoke-v9"
    or plan["plan_sha256"] != "edf243d20e2b8f585bddc0a78fc0b0cee50fa8040206794146dc1e11a9230c95"
    or plan["source_lock_sha256"] != "810da04de9fec741ca626406d9c4a9eac6cb3b3c890d990dc8367a9523fabfe4"
    or plan["resource_totals"]["jobs"] != 175
    or plan["resource_totals"]["jobs_by_profile"] != {
        "cpu_baseline_smoke": 150,
        "gpu_rna_atac_smoke": 25,
    }
    or plan["resource_totals"]["gpu_hours_requested"] != 100.0
    or plan["resource_totals"]["cpu_hours_requested"] != 550.0
    or plan["safety"]["sealed_labels"] is not False
    or plan["safety"]["sealed_features"] is not False
    or plan["safety"]["arrays"] is not False
):
    raise SystemExit("frozen candidate plan differs")
included_paths = plan["source_lock"].get("included_paths")
if source_lock(Path(plan["package_root"]), included_paths=included_paths) != plan["source_lock"]:
    raise SystemExit("frozen candidate scoped source lock drifted")
expected = plan["resource_firewall"]["snapshot_sha256"]
started = time.time()
first = resource_snapshot(repository)
second = resource_snapshot(repository)
if first["snapshot_sha256"] != second["snapshot_sha256"] or first["snapshot_sha256"] != expected:
    raise SystemExit("Resource authorities differ from the frozen candidate")
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
    raise SystemExit("Resource authorities have not met the 30-minute stability gate")
write_json_exclusive(stage / "resource_stability_receipt.json", {
    "schema_version": "masld-bench-cobolt-v9-resource-stability-v1",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v9",
    "candidate": candidate.as_posix(),
    "candidate_artifacts_sha256": candidate_sha,
    "plan_sha256": plan["plan_sha256"],
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
    "candidate_artifacts_sha256": candidate_sha,
    "resource_snapshot_sha256": expected,
    "outcomes_opened": False,
    "status": "pass",
}))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
