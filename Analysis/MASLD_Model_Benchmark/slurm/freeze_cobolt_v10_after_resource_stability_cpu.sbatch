#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-freeze-351
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-freeze-351-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-freeze-351-%j.err

# Re-baseline the Cobolt five-seed campaign separation after GPU job 21100294 was
# stopped at run inclusion by a concurrent Resource-authority adoption.
#
# The Resource must be proven quiet BEFORE the new baseline is frozen, otherwise
# re-freezing mid-edit simply reproduces the failure.  Four independent reads of
# the protected snapshot are required to agree: two probes separated by a wait,
# the freeze's own internal read, and one post-freeze read.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

REPOSITORY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${REPOSITORY}/Analysis/MASLD_Model_Benchmark
CAMPAIGN=${ROOT}/config/campaigns/v1_rna_atac_cobolt_5seed_smoke.toml
PRIOR_CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v9--edf243d20e2b8f58
PRIOR_CANDIDATE_SHA=ac24e1a8da4fde1c35aafa1370659155eded36811099b40f3ee39a4fc2866862
INCIDENT=${ROOT}/config/artifacts/incidents/cobolt_v9_resource_authority_firewall_21100294.json
INCIDENT_SHA=0715cf3a04dc34609a8d6ae662cc3a0140c729115d141863f89e7e575c30860a
TARGET=${ROOT}/executions/cobolt-v10-prefreeze-resource-stability-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the Cobolt v10 pre-freeze stability audit.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/cobolt-v10-freeze-${SLURM_JOB_ID}
cd "${ROOT}"
printf '%s  %s\n' "${PRIOR_CANDIDATE_SHA}" "${PRIOR_CANDIDATE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${INCIDENT_SHA}" "${INCIDENT}" | sha256sum --check --strict
PYTHONPATH=src python -m masld_bench.cli registry validate --config-root config
PYTHONPATH=src python - "${REPOSITORY}" "${ROOT}" "${CAMPAIGN}" "${PRIOR_CANDIDATE}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
import time

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.firewall import resource_snapshot
from masld_bench.hashing import sha256_file
from masld_bench.planner import freeze_campaign, load_frozen_plan

repository, root, campaign, prior, stage = (
    Path(value).resolve(strict=True) for value in sys.argv[1:]
)
CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v10"
PRIOR_CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v9"
MINIMUM_QUIET_SECONDS = 1800.0
PROBE_INTERVAL_SECONDS = 300.0

verify_frozen_tree(prior)
prior_plan = load_frozen_plan(prior)
if prior_plan["campaign"]["campaign_id"] != PRIOR_CAMPAIGN_ID:
    raise SystemExit("superseded candidate identity differs")

# --- stability guard -------------------------------------------------------
started = time.time()
first = resource_snapshot(repository)
newest_mtime = 0.0
newest_path = None
for record in first["files"]:
    observed = (repository / record["path"]).stat(follow_symlinks=False).st_mtime
    if observed > newest_mtime:
        newest_mtime, newest_path = observed, record["path"]
if newest_path is None:
    raise SystemExit("Resource authority snapshot has no protected files")
quiet_seconds = started - newest_mtime
if quiet_seconds < MINIMUM_QUIET_SECONDS:
    raise SystemExit(
        "Resource authorities are still moving: newest protected file "
        f"{newest_path} is only {quiet_seconds:.0f}s old; refusing to re-freeze mid-edit"
    )
time.sleep(PROBE_INTERVAL_SECONDS)
second = resource_snapshot(repository)
if first["snapshot_sha256"] != second["snapshot_sha256"]:
    changed = sorted(
        path
        for path in {item["path"] for item in first["files"]}
        | {item["path"] for item in second["files"]}
        if {item["path"]: item for item in first["files"]}.get(path)
        != {item["path"]: item for item in second["files"]}.get(path)
    )
    raise SystemExit(
        "Resource authorities changed between the two probes; refusing to re-freeze: "
        + ", ".join(changed[:40])
    )
if first["snapshot_sha256"] == prior_plan["resource_firewall"]["snapshot_sha256"]:
    raise SystemExit(
        "current Resource snapshot already equals the superseded v9 baseline; "
        "no re-baseline is warranted and the v9 failure is unexplained"
    )

# --- freeze the revision ---------------------------------------------------
target = freeze_campaign(
    campaign,
    config_root=root / "config",
    package_root=root,
    output_root=root / "candidates",
)
plan = load_frozen_plan(target)
verify_frozen_tree(target)
if plan["campaign"]["campaign_id"] != CAMPAIGN_ID:
    raise SystemExit("frozen campaign identity differs")
if plan["resource_firewall"]["snapshot_sha256"] != first["snapshot_sha256"]:
    raise SystemExit("freeze captured a different Resource snapshot than the verified probes")

third = resource_snapshot(repository)
if third["snapshot_sha256"] != first["snapshot_sha256"]:
    raise SystemExit("Resource authorities changed during the freeze itself")

# --- revision scope: separation baseline only --------------------------------
if plan["source_lock"] != prior_plan["source_lock"]:
    raise SystemExit("revision scope violated: scoped source lock changed")
for key in (
    "runtime_locks",
    "resource_profiles",
    "control_plane_lock",
    "dataset_views",
    "splits",
):
    if key in prior_plan and plan.get(key) != prior_plan.get(key):
        raise SystemExit(f"revision scope violated: {key} changed")
if plan["resource_totals"] != prior_plan["resource_totals"]:
    raise SystemExit("revision scope violated: resource totals changed")
if plan["safety"] != prior_plan["safety"]:
    raise SystemExit("revision scope violated: safety flags changed")

def run_shape(source):
    """Every scientific field of every run, minus the two identity fields that a
    revision is expected to change: the campaign name and the run hash derived
    from it, plus the two digests that necessarily move with the revision text.

    campaign_file_sha256 is the digest of the campaign TOML itself.
    registry_snapshot_sha256 is the aggregate digest of the 14 explicit registry
    files, and the campaign TOML is one of those 14 -- so editing the campaign
    file to declare the revision necessarily changes it too.  Exempting only
    campaign_file_sha256 made this guard self-defeating: the edit a revision
    REQUIRES was rejected as a scope violation.  Both are exempted here, and the
    registry is instead checked file-by-file immediately below, which is
    strictly stronger than comparing the aggregate hash."""
    rows = []
    for run in source["runs"]:
        row = {k: v for k, v in run.items() if k not in ("run_id", "campaign_id")}
        metadata = dict(row["metadata"])
        metadata.pop("campaign_file_sha256", None)
        metadata.pop("registry_snapshot_sha256", None)
        row["metadata"] = metadata
        rows.append(json.dumps(row, sort_keys=True))
    return sorted(rows)

def registry_files(source):
    return {item["path"]: item["sha256"] for item in source["registry_snapshot"]["files"]}

# The registry snapshot may differ ONLY in the campaign file itself.  Every
# other registry member -- datasets, dataset views, splits, tasks, models,
# runtimes, resources, evaluation requirements -- must be byte-identical to the
# superseded revision, or this is not a firewall-baseline-only revision.
CAMPAIGN_REGISTRY_MEMBER = "campaigns/v1_rna_atac_cobolt_5seed_smoke.toml"
new_registry = registry_files(plan)
old_registry = registry_files(prior_plan)
if set(new_registry) != set(old_registry):
    raise SystemExit(
        "revision scope violated: registry membership changed; added="
        + ",".join(sorted(set(new_registry) - set(old_registry)))
        + " removed="
        + ",".join(sorted(set(old_registry) - set(new_registry)))
    )
registry_changed = sorted(
    path for path in new_registry if new_registry[path] != old_registry[path]
)
if registry_changed != [CAMPAIGN_REGISTRY_MEMBER]:
    raise SystemExit(
        "revision scope violated: registry files other than the campaign file "
        "changed: " + ", ".join(p for p in registry_changed if p != CAMPAIGN_REGISTRY_MEMBER)
    )

if run_shape(plan) != run_shape(prior_plan):
    raise SystemExit(
        "revision scope violated: a seed, fold, model, hyperparameter, dataset "
        "binding, split, runtime, resource profile, evaluator, or endpoint differs "
        "from the superseded revision"
    )
digests = {run["metadata"]["campaign_file_sha256"] for run in plan["runs"]}
prior_digests = {run["metadata"]["campaign_file_sha256"] for run in prior_plan["runs"]}
if len(digests) != 1 or digests == prior_digests:
    raise SystemExit("revision campaign-file digest is not a single new value")
prior_ids = {str(run["run_id"]) for run in prior_plan["runs"]}
new_ids = {str(run["run_id"]) for run in plan["runs"]}
if prior_ids & new_ids:
    raise SystemExit("revision reuses superseded run identities")
if plan["plan_sha256"] == prior_plan["plan_sha256"]:
    raise SystemExit("revision is not distinct from the superseded plan")

cobolt = sorted(
    (
        (int(run["seed"]), int(run["fold"]))
        for run in plan["runs"]
        if run.get("model_id") == "cobolt"
        and run.get("resource_profile") == "gpu_rna_atac_smoke"
    )
)
expected = [(seed, fold) for seed in (1103, 2107, 3109, 4111, 5113) for fold in range(5)]
if cobolt != expected:
    raise SystemExit("frozen v10 does not carry the unchanged five-seed by five-fold grid")

receipt = {
    "schema_version": "masld-bench-cobolt-v10-prefreeze-resource-stability-v1",
    "status": "pass",
    "campaign_id": CAMPAIGN_ID,
    "supersedes_campaign_id": PRIOR_CAMPAIGN_ID,
    "supersedes_failed_job_ids": [21100294],
    "candidate": target.as_posix(),
    "candidate_artifacts_sha256": sha256_file(target / "ARTIFACTS.json"),
    "plan_sha256": plan["plan_sha256"],
    "superseded_candidate": prior.as_posix(),
    "superseded_candidate_artifacts_sha256": sha256_file(prior / "ARTIFACTS.json"),
    "superseded_plan_sha256": prior_plan["plan_sha256"],
    "superseded_resource_snapshot_sha256": prior_plan["resource_firewall"]["snapshot_sha256"],
    "resource_snapshot_sha256": first["snapshot_sha256"],
    "authority_root_count": len(first["authority_roots"]),
    "authority_file_count": len(first["files"]),
    "newest_protected_file": newest_path,
    "newest_protected_file_mtime_unix": newest_mtime,
    "observed_quiet_seconds_at_start": quiet_seconds,
    "minimum_quiet_seconds": MINIMUM_QUIET_SECONDS,
    "probe_interval_seconds": PROBE_INTERVAL_SECONDS,
    "agreeing_snapshot_reads": 4,
    "read_order": ["probe_before_wait", "probe_after_wait", "freeze_internal", "probe_after_freeze"],
    "source_lock_unchanged": True,
    "resource_totals_unchanged": True,
    "safety_flags_unchanged": True,
    "run_grid_unchanged": True,
    "run_identities_disjoint_from_superseded": True,
    "firewall_watched_roots_reduced": False,
    "outcomes_opened": False,
    "metrics_calculated": False,
    "gpu_jobs_submitted": 0,
}
write_json_exclusive(stage / "resource_stability_receipt.json", receipt)
print(
    freeze_tree(
        stage,
        {
            "artifact_class": "cobolt_v10_prefreeze_resource_stability",
            "campaign_id": CAMPAIGN_ID,
            "candidate_artifacts_sha256": receipt["candidate_artifacts_sha256"],
            "resource_snapshot_sha256": receipt["resource_snapshot_sha256"],
            "outcomes_opened": False,
            "status": "pass",
        },
    )
)
print(json.dumps({"candidate": target.as_posix(), "plan_sha256": plan["plan_sha256"]}, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
