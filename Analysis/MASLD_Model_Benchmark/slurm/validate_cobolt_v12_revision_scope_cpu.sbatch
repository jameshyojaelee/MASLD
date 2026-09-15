#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-357
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-check-357-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-check-357-%j.err

# Revision-scope validation of the ALREADY FROZEN Cobolt v12 candidate.
#
# Freeze job 21110827 wrote candidates/v1-rna-atac-cobolt-5seed-smoke-v12--3972f43112b8ecc4
# with COMPLETE, then FAILED in its post-freeze scope check on a variable-name bug
# of mine: the block subscripted `prior` (a PosixPath) instead of `prior_plan`.
# Nine sites, not one. Validating the existing candidate instead of re-freezing
# saves a third 34-minute walk of 153,568 files and, more importantly, removes a
# third 34-minute window in which a Fig4 adoption could land.
#
# STANDING RULE THIS COST US: any logic that runs AFTER an expensive step must be
# exercised against existing frozen output files BEFORE the expensive step is paid
# for. Every assertion below operates on two already-frozen plan dicts, so the
# whole block was dry-runnable in seconds against v10/v11/v12 on disk. Doing that
# found THREE bugs, not one:
#   1. seven `prior[...]` subscripts (the reported TypeError)
#   2. two bare `prior` uses that a `prior[` grep missed - caught only because a
#      word-boundary regex was used after the dry-run still failed
#   3. a v11-wrapper remnant comparing the prior GPU grid on OLD_PROFILE, when a
#      separation-only revision leaves BOTH v11 and v12 on the B6K profile
# Bugs 2 and 3 would each have cost their own 34-minute freeze.
#
# Dry-run verdicts before this job was written:
#   v12 vs v11 -> PASS          (target)
#   v11 vs v10 -> REJECTED      (runtime revision; negative control)
#   v12 vs v12 -> REJECTED      (snapshot unchanged; negative control)

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

REPOSITORY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${REPOSITORY}/Analysis/MASLD_Model_Benchmark
V12=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v12--3972f43112b8ecc4
V11=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v11--e17be3dab9ba48cb
RUNTIME_TREE=${ROOT}/executions/environments/cobolt-native-torch2.8.0-cu128-21107432
RUNTIME_TREE_SHA=e7e56c0ee207843b78dd9bbb1e3619fc965fc7836a4d3b0131c58b6f4e58049c
VALIDATION=${ROOT}/executions/cobolt-b6k-runtime-validation-21107495
VALIDATION_SHA=9ab79607dc34c059485ed04e4c2fef005713eddd8003cf4d234a5ecdf16afa7b
TARGET=${ROOT}/executions/cobolt-v12-revision-scope-validation-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the Cobolt v12 scope validation.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/cobolt-v12-scope-${SLURM_JOB_ID}
cd "${ROOT}"
sha256sum --check --strict <<EOF
${RUNTIME_TREE_SHA}  ${RUNTIME_TREE}/ARTIFACTS.json
${VALIDATION_SHA}  ${VALIDATION}/ARTIFACTS.json
EOF
( cd "${RUNTIME_TREE}" && sha256sum --check --strict SHA256SUMS )

PYTHONPATH=src python - "${REPOSITORY}" "${V12}" "${V11}" "${STAGE}" "${RUNTIME_TREE}" "${VALIDATION}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.firewall import resource_snapshot
from masld_bench.hashing import sha256_file
from masld_bench.planner import load_frozen_plan

repository, v12p, v11p, stage, runtime_tree, validation = (
    Path(v).resolve(strict=True) for v in sys.argv[1:]
)
NEW_RUNTIME_ID = "gpu_rna_atac_torch_b6k"
OLD_RUNTIME_ID = "gpu_rna_atac_torch_smoke"
NEW_PROFILE, OLD_PROFILE = "gpu_rna_atac_b6k", "gpu_rna_atac_smoke"
root = v12p.parent.parent

for tree in (v12p, v11p, validation):
    verify_frozen_tree(tree)
plan, prior_plan = load_frozen_plan(v12p), load_frozen_plan(v11p)
if plan["campaign"]["campaign_id"] != "v1-rna-atac-cobolt-5seed-smoke-v12":
    raise SystemExit("v12 campaign id differs")
if prior_plan["campaign"]["campaign_id"] != "v1-rna-atac-cobolt-5seed-smoke-v11":
    raise SystemExit("v11 campaign id differs")

vr = json.loads((validation / "validation.json").read_text())
for c in ("forward_backward_passed", "query_batch_invariance_passed",
          "tensor_state_restore_passed", "cross_process_determinism_passed",
          "environment_contract_passed"):
    if vr.get(c) is not True:
        raise SystemExit(f"B6K runtime validation check failed: {c}")
if vr.get("l40s_validation_reused") is not False:
    raise SystemExit("B6K validation reused the L40S validation")

try:
    # --- v12 is FIREWALL-BASELINE-ONLY: the opposite scope from v11 --------------
    # v11 was required to change the fixed source record, fixed runtime and profile. v12 must
    # change NOTHING except the protected Resource snapshot, so every one of those is
    # asserted UNCHANGED here.
    if plan["source_lock"] != prior_plan["source_lock"]:
        raise SystemExit("scope violated: v12 must not change the scoped source lock")
    if plan["runtime_locks"] != prior_plan["runtime_locks"]:
        raise SystemExit("scope violated: v12 must not change any runtime lock")
    if plan["resource_profiles"] != prior_plan["resource_profiles"]:
        raise SystemExit("scope violated: v12 must not change any resource profile")
    for key in ("dataset_locks", "safety", "resource_totals"):
        if plan.get(key) != prior_plan.get(key):
            raise SystemExit(f"scope violated: {key} changed")
    if NEW_RUNTIME_ID not in plan["runtime_locks"] or OLD_RUNTIME_ID in plan["runtime_locks"]:
        raise SystemExit("v12 does not bind the B6K runtime exclusively")
    if plan["resource_profiles"][NEW_PROFILE]["accelerator"] != "b6k":
        raise SystemExit("v12 GPU profile does not request the B6K accelerator")
    if str(runtime_tree) not in " ".join(str(x) for x in plan["runtime_locks"][NEW_RUNTIME_ID]["command_prefix"]):
        raise SystemExit("v12 runtime lock does not point at the validated B6K runtime tree")

    # The one thing that MUST move: the protected Resource baseline.
    if plan["resource_firewall"]["snapshot_sha256"] == prior_plan["resource_firewall"]["snapshot_sha256"]:
        raise SystemExit(
            "v12 captured the same Resource snapshot as v11; no re-baseline was warranted "
            "and the v11 firewall abort is unexplained"
        )
    if plan["plan_sha256"] == prior_plan["plan_sha256"]:
        raise SystemExit("revision is not distinct from the superseded plan")
    if set(str(r["run_id"]) for r in plan["runs"]) & set(str(r["run_id"]) for r in prior_plan["runs"]):
        raise SystemExit("revision reuses superseded run identities")

    # --- everything scientific, INCLUDING the environment inputs, must be identical
    # v11->v12 changes no runtime, so unlike v11 there must be ZERO environment-role
    # deltas. The environment input is compared here rather than set aside.
    _EXEMPT_RUN_FIELDS = {"run_id", "campaign_id"}
    _SCI = ("task_id", "model_id", "dataset_ids", "seed", "fold", "split_id", "stage",
            "action", "adaptation_regime", "hyperparameters", "checkpoint", "inputs",
            "immutable_inputs", "adapter_command", "runtime_id", "resource_profile",
            "runtime_registry_sha256", "code_lock_sha256")

    def ident(run):
        return (run.get("model_id"), run.get("seed"), run.get("fold"),
                run.get("stage"), tuple(run.get("action") or []))

    new_runs = {ident(r): r for r in plan["runs"]}
    old_runs = {ident(r): r for r in prior_plan["runs"]}
    if set(new_runs) != set(old_runs):
        raise SystemExit("scope violated: the set of scientific runs changed")
    if len(plan["runs"]) != 175:
        raise SystemExit(f"expected 175 runs, found {len(plan['runs'])}")
    for k in sorted(set(new_runs) & set(old_runs), key=str):
        a, b = old_runs[k], new_runs[k]
        for f in _SCI:
            if a.get(f) != b.get(f):
                raise SystemExit(f"scope violated: field {f} changed for run {k}")
        am = dict(a.get("metadata") or {})
        bm = dict(b.get("metadata") or {})
        # campaign_file_sha256 and registry_snapshot_sha256 necessarily move with the
        # revision text; the registry is checked file-by-file below instead.
        for drop in ("campaign_file_sha256", "registry_snapshot_sha256"):
            am.pop(drop, None); bm.pop(drop, None)
        if am != bm:
            differing = sorted(x for x in set(am) | set(bm) if am.get(x) != bm.get(x))
            raise SystemExit(f"scope violated: metadata {differing} changed for run {k}")

    # --------------------------------------------------------------------- registry
    def registry_files(source):
        return {i["path"]: i["sha256"] for i in source["registry_snapshot"]["files"]}

    new_registry, old_registry = registry_files(plan), registry_files(prior_plan)
    if set(new_registry) != set(old_registry):
        raise SystemExit("scope violated: registry membership changed")
    changed = sorted(p for p in new_registry if new_registry[p] != old_registry[p])
    if changed != ["campaigns/v1_rna_atac_cobolt_5seed_smoke.toml"]:
        raise SystemExit(f"scope violated: registry files beyond the campaign changed: {changed}")

    models_toml = (root / "config" / "models" / "multimodal_integration.toml").read_text()
    if 'runtime_id = "gpu_rna_atac_torch_smoke"' not in models_toml:
        raise SystemExit("scope violated: scpair no longer binds the L40S runtime")
    if models_toml.count('runtime_id = "gpu_rna_atac_torch_b6k"') != 1:
        raise SystemExit("scope violated: more than one model binds B6K")

    # --- the unchanged five-seed by five-fold grid, on the new profile ----------
    grid = sorted(
        (int(run["seed"]), int(run["fold"]))
        for run in plan["runs"]
        if run.get("model_id") == "cobolt" and run.get("resource_profile") == NEW_PROFILE
    )
    expected_grid = [(s, f) for s in (1103, 2107, 3109, 4111, 5113) for f in range(5)]
    if grid != expected_grid:
        raise SystemExit("frozen v12 does not carry the unchanged five-seed by five-fold grid")
    prior_grid = sorted(
        (int(run["seed"]), int(run["fold"]))
        for run in prior_plan["runs"]
        # v11 and v12 are BOTH on the B6K profile: a separation-only revision does not
        # move the profile. Comparing against OLD_PROFILE here was a v11-wrapper remnant.
        if run.get("model_id") == "cobolt" and run.get("resource_profile") == NEW_PROFILE
    )
    if grid != prior_grid:
        raise SystemExit("v12 GPU grid differs from v11")
    cpu_runs = [r for r in plan["runs"] if r.get("resource_profile") != NEW_PROFILE]
    if len(cpu_runs) != 150:
        raise SystemExit(f"expected 150 CPU control runs, found {len(cpu_runs)}")


except SystemExit as error:
    raise SystemExit(f"v12 revision scope FAILED: {error}")

observed = resource_snapshot(repository)
frozen = plan["resource_firewall"]["snapshot_sha256"]
still_current = observed["snapshot_sha256"] == frozen

receipt = {
    "schema_version": "masld-bench-cobolt-v12-revision-scope-validation-v1",
    "status": "passed_firewall_baseline_only_revision",
    "candidate": v12p.as_posix(),
    "candidate_manifest_sha256": sha256_file(v12p / "ARTIFACTS.json"),
    "plan_sha256": plan["plan_sha256"],
    "supersedes_candidate": v11p.as_posix(),
    "supersedes_plan_sha256": prior_plan["plan_sha256"],
    "revision_class": "firewall_baseline_only",
    "runs_total": len(plan["runs"]),
    "gpu_runs_on_b6k": 25,
    "cpu_control_runs": 150,
    "scientific_content_identical_to_v11": True,
    "resource_firewall_snapshot_at_freeze": frozen,
    "resource_firewall_snapshot_now": observed["snapshot_sha256"],
    "resource_firewall_still_current": still_current,
    "prior_firewall_snapshot_sha256": prior_plan["resource_firewall"]["snapshot_sha256"],
    "why_v11_was_superseded": (
        "v11 had NO defect of its own - its scope validation passed as "
        "passed_execution_runtime_only_revision. It was invalidated because a protected "
        "Resource authority moved in the ~10-minute window after job 21109214 recorded "
        "resource_firewall_unchanged_since_freeze: true. figures/main/fig4_singlecell_"
        "programs/captions.md was rewritten at 12:53:11. v11 remains a valid frozen record."
    ),
    "freeze_job_failed_after_writing_the_candidate": (
        "Job 21110827 wrote this candidate with COMPLETE and then failed in its own "
        "post-freeze scope check on a variable-name bug. Third consecutive freeze where "
        "the candidate landed and the job then failed, so the reliable signal is COMPLETE "
        "on the tree plus a terminal job state - NEVER the candidate directory alone."
    ),
    "dry_run_before_submission": {
        "v12_vs_v11": "PASS", "v11_vs_v10": "REJECTED", "v12_vs_v12": "REJECTED",
        "bugs_found_by_dry_run": 3,
    },
    "runtime_tree_artifacts_sha256": sha256_file(runtime_tree / "ARTIFACTS.json"),
    "b6k_validation_artifacts_sha256": sha256_file(validation / "ARTIFACTS.json"),
    "l40s_validation_reused": False,
}
if not still_current:
    receipt["status"] = "failed_firewall_moved_since_freeze"
(stage / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
freeze_tree(stage, {"artifact_class": "cobolt_v12_revision_scope_validation"})
print(json.dumps(receipt, indent=2, sort_keys=True))
if not still_current:
    raise SystemExit("Resource firewall moved since the v12 freeze; do not launch")
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Cobolt v12 revision scope validated: %s\n' "${TARGET}"
