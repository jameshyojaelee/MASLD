#!/usr/bin/env bash
#SBATCH --job-name=model-check-352
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Validate that the ALREADY-FROZEN Cobolt v10 candidate is a firewall-baseline-only
# revision of v9.
#
# Job 21104557 froze the v10 candidate successfully and then failed its own scope
# check.  Diagnosis: the check compared metadata.registry_snapshot_sha256, but the
# campaign TOML is one of the 14 explicit registry files, so the edit that DECLARES
# the revision necessarily changes that aggregate digest.  The guard rejected the
# very edit a revision requires.  This validator exempts that aggregate digest and
# instead compares the registry file-by-file, which is strictly stronger.
#
# This validator does NOT re-freeze, does NOT publish to the GPU queue, and does
# NOT submit any GPU work.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

REPOSITORY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${REPOSITORY}/Analysis/MASLD_Model_Benchmark
PRIOR=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v9--edf243d20e2b8f58
PRIOR_SHA=ac24e1a8da4fde1c35aafa1370659155eded36811099b40f3ee39a4fc2866862
NEW=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v10--6cbef395fdccc389
TARGET=${ROOT}/executions/cobolt-v10-revision-scope-validation-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the Cobolt v10 scope validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/cobolt-v10-scope-${SLURM_JOB_ID}
cd "${ROOT}"
printf '%s  %s\n' "${PRIOR_SHA}" "${PRIOR}/ARTIFACTS.json" | sha256sum --check --strict

python3 -m pip freeze --all > "${STAGE}/pip_freeze.txt" 2>/dev/null || true
sha256sum \
  slurm/validate_cobolt_v10_revision_scope_cpu.sbatch \
  "${PRIOR}/ARTIFACTS.json" \
  "${NEW}/ARTIFACTS.json" > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src:${ROOT}" python - "$STAGE" "$PRIOR" "$NEW" "$REPOSITORY" <<'PY'
import json, sys
from pathlib import Path
from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.firewall import resource_snapshot

stage, prior_path, new_path, repository = (Path(a) for a in sys.argv[1:5])

verify_frozen_tree(prior_path)
verify_frozen_tree(new_path)
prior_plan = json.loads((prior_path / "plan.json").read_text())
plan = json.loads((new_path / "plan.json").read_text())

findings = {}

def run_shape(source):
    rows = []
    for run in source["runs"]:
        row = {k: v for k, v in run.items() if k not in ("run_id", "campaign_id")}
        metadata = dict(row["metadata"])
        metadata.pop("campaign_file_sha256", None)
        metadata.pop("registry_snapshot_sha256", None)
        row["metadata"] = metadata
        rows.append(json.dumps(row, sort_keys=True))
    return sorted(rows)

findings["run_count_prior"] = len(prior_plan["runs"])
findings["run_count_new"] = len(plan["runs"])
findings["scientific_run_shape_identical"] = run_shape(plan) == run_shape(prior_plan)

def registry_files(source):
    return {i["path"]: i["sha256"] for i in source["registry_snapshot"]["files"]}

CAMPAIGN_MEMBER = "campaigns/v1_rna_atac_cobolt_5seed_smoke.toml"
new_reg, old_reg = registry_files(plan), registry_files(prior_plan)
findings["registry_membership_identical"] = set(new_reg) == set(old_reg)
findings["registry_member_count"] = len(new_reg)
changed = sorted(p for p in new_reg if p in old_reg and new_reg[p] != old_reg[p])
findings["registry_files_changed"] = changed
findings["only_campaign_file_changed"] = changed == [CAMPAIGN_MEMBER]

for key in ("source_lock", "runtime_locks", "resource_profiles",
            "control_plane_lock", "dataset_views", "splits",
            "resource_totals", "safety"):
    if key in prior_plan or key in plan:
        findings[f"identical__{key}"] = plan.get(key) == prior_plan.get(key)

digests = {r["metadata"]["campaign_file_sha256"] for r in plan["runs"]}
prior_digests = {r["metadata"]["campaign_file_sha256"] for r in prior_plan["runs"]}
findings["campaign_digest_single_new_value"] = len(digests) == 1 and digests != prior_digests
findings["no_reused_run_identities"] = not (
    {str(r["run_id"]) for r in prior_plan["runs"]} & {str(r["run_id"]) for r in plan["runs"]}
)
findings["plan_distinct"] = plan["plan_sha256"] != prior_plan["plan_sha256"]
findings["firewall_baseline_changed"] = (
    plan["resource_firewall"]["snapshot_sha256"]
    != prior_plan["resource_firewall"]["snapshot_sha256"]
)

roots = [str(p) for p in prior_plan["resource_firewall"]["authority_roots"]]
observed = resource_snapshot(repository, roots)
findings["frozen_baseline_still_current"] = (
    observed["snapshot_sha256"] == plan["resource_firewall"]["snapshot_sha256"]
)

required = [
    "scientific_run_shape_identical", "registry_membership_identical",
    "only_campaign_file_changed", "campaign_digest_single_new_value",
    "no_reused_run_identities", "plan_distinct", "firewall_baseline_changed",
] + [k for k in findings if k.startswith("identical__")]
failed = [k for k in required if findings.get(k) is not True]

receipt = {
    "schema_version": "masld-bench-cobolt-v10-revision-scope-validation-v1",
    "prior_candidate": str(prior_path),
    "new_candidate": str(new_path),
    "superseded_failed_job_ids": [21100294, 21104557],
    "diagnosis": (
        "Job 21104557 froze the v10 candidate successfully and then failed its own "
        "scope check on metadata.registry_snapshot_sha256. The campaign TOML is one "
        "of the 14 explicit registry files, so editing it to declare the revision "
        "necessarily changes that aggregate digest. The guard therefore rejected the "
        "exact edit a revision requires. Corrected by exempting the aggregate digest "
        "and comparing the registry file-by-file instead, which is strictly stronger."
    ),
    "findings": findings,
    "failed_checks": failed,
    "status": "passed_firewall_baseline_only_revision" if not failed else "failed_scope_violation",
    "gpu_submitted": False,
    "queue_published": False,
    "refroze_candidate": False,
    "outcomes_read": False,
    "metrics_calculated": False,
}
write_json_exclusive(stage / "scope_validation.json", receipt)
freeze_tree(stage, {
    "artifact_class": "cobolt_v10_revision_scope_validation",
    "status": receipt["status"],
    "outcomes_read": False,
    "metrics_calculated": False,
    "gpu_submitted": False,
    "queue_published": False,
})
print(json.dumps({"status": receipt["status"], "failed_checks": failed}, sort_keys=True))
if failed:
    raise SystemExit("v10 revision scope validation failed: " + ", ".join(failed))
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
