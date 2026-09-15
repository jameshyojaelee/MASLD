#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-freeze-352
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-freeze-352-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-freeze-352-%j.err

# Freeze Cobolt campaign v11: an EXECUTION-RUNTIME revision that moves the 25 GPU
# fits from the saturated L40S queue to the idle B6K pool.
#
# This is deliberately NOT a firewall-baseline-only revision, so the v10 scope
# guard is inverted here.  v11 is REQUIRED to change the scoped fixed source record, the
# GPU fixed runtime, the GPU resource profile, and the registry membership; it is
# FORBIDDEN to change any scientific quantity.  Both directions are asserted
# below, so a revision that silently altered a seed, fold, hyperparameter,
# dataset binding, split, or evaluator would fail here rather than at analysis.
#
# The Resource stability guard is unchanged from the v10 freeze, because that is
# what job 21100294 needed and did not have: all five authority documents were
# edited earlier today, and a concurrent adoption aborts every run of a bundle.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

REPOSITORY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${REPOSITORY}/Analysis/MASLD_Model_Benchmark
CAMPAIGN=${ROOT}/config/campaigns/v1_rna_atac_cobolt_5seed_smoke.toml
PRIOR_CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v10--6cbef395fdccc389
RUNTIME_TREE=${ROOT}/executions/environments/cobolt-native-torch2.8.0-cu128-21107432
RUNTIME_TREE_SHA=e7e56c0ee207843b78dd9bbb1e3619fc965fc7836a4d3b0131c58b6f4e58049c
VALIDATION=${ROOT}/executions/cobolt-b6k-runtime-validation-21107495
VALIDATION_SHA=9ab79607dc34c059485ed04e4c2fef005713eddd8003cf4d234a5ecdf16afa7b
ADAPTER_TEST=${ROOT}/executions/cobolt-b6k-adapter-test-21107681
TARGET=${ROOT}/executions/cobolt-v11-prefreeze-resource-stability-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the Cobolt v11 pre-freeze stability audit.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/cobolt-v11-freeze-${SLURM_JOB_ID}
cd "${ROOT}"

# The runtime and its two on-hardware validations must be intact before the
# campaign that binds them is frozen.
sha256sum --check --strict <<EOF
${RUNTIME_TREE_SHA}  ${RUNTIME_TREE}/ARTIFACTS.json
${VALIDATION_SHA}  ${VALIDATION}/ARTIFACTS.json
EOF
# The runtime tree carries a conda prefix, so it is pinned by its own SHA256SUMS
# manifest rather than by freeze_tree(), whose inventory_tree() rejects symlinks.
( cd "${RUNTIME_TREE}" && sha256sum --check --strict SHA256SUMS )

PYTHONPATH=src python -m masld_bench.cli registry validate --config-root config

PYTHONPATH=src python - \
  "${REPOSITORY}" "${ROOT}" "${CAMPAIGN}" "${PRIOR_CANDIDATE}" "${STAGE}" \
  "${RUNTIME_TREE}" "${VALIDATION}" "${ADAPTER_TEST}" <<'PY'
from pathlib import Path
import json
import sys
import time

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.firewall import resource_snapshot
from masld_bench.hashing import sha256_file
from masld_bench.planner import freeze_campaign, load_frozen_plan

(
    repository, root, campaign, prior, stage,
    runtime_tree, validation, adapter_test,
) = (Path(value).resolve(strict=True) for value in sys.argv[1:])

CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v11"
PRIOR_CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v10"
NEW_RUNTIME_ID = "gpu_rna_atac_torch_b6k"
OLD_RUNTIME_ID = "gpu_rna_atac_torch_smoke"
NEW_PROFILE = "gpu_rna_atac_b6k"
OLD_PROFILE = "gpu_rna_atac_smoke"
MINIMUM_QUIET_SECONDS = 1800.0
PROBE_INTERVAL_SECONDS = 300.0

verify_frozen_tree(validation)
if not (runtime_tree / "SHA256SUMS").is_file():
    raise SystemExit("runtime tree has no SHA256SUMS manifest")
adapter_result = json.loads((adapter_test / "result.json").read_text())
validation_receipt = json.loads((validation / "validation.json").read_text())
if validation_receipt.get("status") != "pass":
    raise SystemExit("B6K runtime validation did not pass")
if validation_receipt.get("l40s_validation_reused") is not False:
    raise SystemExit("B6K validation must not reuse the L40S validation")
for check in (
    "forward_backward_passed",
    "query_batch_invariance_passed",
    "tensor_state_restore_passed",
    "cross_process_determinism_passed",
    "environment_contract_passed",
):
    if validation_receipt.get(check) is not True:
        raise SystemExit(f"B6K validation check did not pass: {check}")

verify_frozen_tree(prior)
prior_plan = load_frozen_plan(prior)
if prior_plan["campaign"]["campaign_id"] != PRIOR_CAMPAIGN_ID:
    raise SystemExit("superseded candidate identity differs")

# --- Resource stability guard (unchanged from the v10 freeze) ---------------
started = time.time()
first = resource_snapshot(repository)
newest_mtime, newest_path = 0.0, None
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
        f"{newest_path} is only {quiet_seconds:.0f}s old; refusing to freeze mid-edit"
    )
time.sleep(PROBE_INTERVAL_SECONDS)
second = resource_snapshot(repository)
if first["snapshot_sha256"] != second["snapshot_sha256"]:
    before = {item["path"]: item for item in first["files"]}
    after = {item["path"]: item for item in second["files"]}
    changed = sorted(p for p in set(before) | set(after) if before.get(p) != after.get(p))
    raise SystemExit(
        "Resource authorities changed between the two probes; refusing to freeze: "
        + ", ".join(changed[:40])
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
    raise SystemExit("freeze captured a different Resource snapshot than the probes")
third = resource_snapshot(repository)
if third["snapshot_sha256"] != first["snapshot_sha256"]:
    raise SystemExit("Resource authorities changed during the freeze itself")

# --- what a runtime revision MUST change -----------------------------------
if plan["source_lock"] == prior_plan["source_lock"]:
    raise SystemExit("v11 must change the scoped source lock; the adapter guard was edited")
if plan["plan_sha256"] == prior_plan["plan_sha256"]:
    raise SystemExit("revision is not distinct from the superseded plan")
if set(str(r["run_id"]) for r in plan["runs"]) & set(
    str(r["run_id"]) for r in prior_plan["runs"]
):
    raise SystemExit("revision reuses superseded run identities")
if NEW_RUNTIME_ID not in plan["runtime_locks"] or OLD_RUNTIME_ID in plan["runtime_locks"]:
    raise SystemExit("v11 must bind the B6K runtime and drop the L40S runtime")
if plan["resource_profiles"][NEW_PROFILE]["accelerator"] != "b6k":
    raise SystemExit("v11 GPU profile does not request the B6K accelerator")
lock = plan["runtime_locks"][NEW_RUNTIME_ID]
if str(runtime_tree) not in " ".join(str(x) for x in lock["command_prefix"]):
    raise SystemExit("v11 runtime lock does not point at the validated B6K runtime tree")

# --- what a runtime revision MUST NOT change -------------------------------
# Compare the SCIENTIFIC fields of every run explicitly rather than diffing the
# whole run dict.  A runtime revision necessarily moves the execution bindings
# (interpreter path, environment pin output file, model-registry digests, code selection record,
# run id), and a whole-dict diff would reject the very change v11 exists to make.
# Each moved binding is asserted separately: the fixed runtime above, the registry
# file-by-file check and scpair guard below.
_SCIENCE_FIELDS = (
    "task_id", "model_id", "dataset_ids", "seed", "fold", "split_id",
    "stage", "action", "adaptation_regime", "hyperparameters", "checkpoint",
)
_SCIENCE_METADATA = (
    "adapter_request", "adapter_timeout_seconds", "dataset_routing",
    "dataset_routing_wave", "development_prediction_dataset_ids",
    "evaluator_contract_sha256", "evaluator_parameters", "fit_dataset_ids",
    "input_owner_by_role", "open_champion_eligible",
    "prediction_first_stress_dataset_ids", "primary_endpoint_scoring_allowed",
    "primary_evaluator_id", "sealed_prediction_dataset_ids", "split_contract",
    "split_contract_sha256", "split_registry_sha256", "task_contract_sha256",
    "task_registry_sha256", "wave",
)

def science_shape(source):
    rows = []
    for run in source["runs"]:
        metadata = run.get("metadata") or {}
        row = {k: run.get(k) for k in _SCIENCE_FIELDS}
        row["metadata"] = {k: metadata.get(k) for k in _SCIENCE_METADATA}
        # Data inputs must be identical; only the environment output file may move.
        row["inputs"] = sorted(
            json.dumps(i, sort_keys=True)
            for i in (run.get("inputs") or [])
            if not str(i.get("role", "")).startswith("environment:")
        )
        rows.append(json.dumps(row, sort_keys=True))
    return sorted(rows)

if science_shape(plan) != science_shape(prior_plan):
    raise SystemExit(
        "scope violated: a seed, fold, model, hyperparameter, dataset binding, "
        "split, evaluator, endpoint, adapter timeout, or data input differs from v10"
    )
# Exactly one environment output file role may have moved, and only to the B6K selection record.
new_env_roles = {
    str(i.get("role"))
    for run in plan["runs"] for i in (run.get("inputs") or [])
    if str(i.get("role", "")).startswith("environment:")
}
if new_env_roles != {f"environment:{NEW_RUNTIME_ID}", "environment:cpu_baseline_smoke"}:
    raise SystemExit(f"unexpected environment artifact roles: {sorted(new_env_roles)}")
for key in ("dataset_locks", "safety"):
    if prior_plan.get(key) != plan.get(key):
        raise SystemExit(f"scope violated: {key} changed")

def registry_files(source):
    return {i["path"]: i["sha256"] for i in source["registry_snapshot"]["files"]}

new_registry, old_registry = registry_files(plan), registry_files(prior_plan)
expected_added = {f"runtimes/{NEW_RUNTIME_ID}.toml"}
expected_removed = {f"runtimes/{OLD_RUNTIME_ID}.toml"}
if set(new_registry) - set(old_registry) != expected_added:
    raise SystemExit("scope violated: unexpected registry additions")
if set(old_registry) - set(new_registry) != expected_removed:
    raise SystemExit("scope violated: unexpected registry removals")
# models/multimodal_integration.toml carries cobolt's runtime_id and
# environment_artifact role, which the planner cross-checks against the campaign
# runtime.  It also carries scpair, which must stay on the L40S runtime; that is
# asserted separately below.
ALLOWED_CHANGED = {
    "campaigns/v1_rna_atac_cobolt_5seed_smoke.toml",
    "models/multimodal_integration.toml",
    "resources.toml",
}
changed = {
    p for p in set(new_registry) & set(old_registry)
    if new_registry[p] != old_registry[p]
}
if changed - ALLOWED_CHANGED:
    raise SystemExit(
        "scope violated: registry files beyond the campaign and resources "
        "changed: " + ", ".join(sorted(changed - ALLOWED_CHANGED))
    )

# scpair must not have been dragged onto B6K by the cobolt edit.
models_toml = (root / "config" / "models" / "multimodal_integration.toml").read_text()
if 'runtime_id = "gpu_rna_atac_torch_smoke"' not in models_toml:
    raise SystemExit("scope violated: scpair no longer binds the L40S runtime")
if models_toml.count('runtime_id = "gpu_rna_atac_torch_b6k"') != 1:
    raise SystemExit("scope violated: more than one model was moved to B6K")

# --- the unchanged five-seed by five-fold grid, on the new profile ----------
grid = sorted(
    (int(run["seed"]), int(run["fold"]))
    for run in plan["runs"]
    if run.get("model_id") == "cobolt" and run.get("resource_profile") == NEW_PROFILE
)
expected_grid = [(s, f) for s in (1103, 2107, 3109, 4111, 5113) for f in range(5)]
if grid != expected_grid:
    raise SystemExit("frozen v11 does not carry the unchanged five-seed by five-fold grid")
prior_grid = sorted(
    (int(run["seed"]), int(run["fold"]))
    for run in prior_plan["runs"]
    if run.get("model_id") == "cobolt" and run.get("resource_profile") == OLD_PROFILE
)
if grid != prior_grid:
    raise SystemExit("v11 GPU grid differs from v10")
cpu_runs = [r for r in plan["runs"] if r.get("resource_profile") != NEW_PROFILE]
if len(cpu_runs) != 150:
    raise SystemExit(f"expected 150 CPU control runs, found {len(cpu_runs)}")

receipt = {
    "schema_version": "masld-bench-cobolt-v11-prefreeze-resource-stability-v1",
    "status": "pass",
    "campaign_id": CAMPAIGN_ID,
    "supersedes_campaign_id": PRIOR_CAMPAIGN_ID,
    "revision_class": "execution_runtime_only",
    "candidate": target.as_posix(),
    "candidate_manifest_sha256": sha256_file(target / "ARTIFACTS.json"),
    "plan_sha256": plan["plan_sha256"],
    "prior_plan_sha256": prior_plan["plan_sha256"],
    "runtime_id": NEW_RUNTIME_ID,
    "resource_profile": NEW_PROFILE,
    "accelerator": "b6k",
    "runtime_tree": runtime_tree.as_posix(),
    "runtime_tree_artifacts_sha256": sha256_file(runtime_tree / "ARTIFACTS.json"),
    "runtime_tree_sha256sums_sha256": sha256_file(runtime_tree / "SHA256SUMS"),
    "b6k_validation": validation.as_posix(),
    "b6k_validation_artifacts_sha256": sha256_file(validation / "ARTIFACTS.json"),
    "b6k_adapter_test": adapter_test.as_posix(),
    "b6k_adapter_test_result_sha256": sha256_file(adapter_test / "result.json"),
    "b6k_adapter_test_equation_parity_passed": adapter_result.get("equation_parity_passed"),
    "b6k_adapter_test_held_atac_exposed": adapter_result.get("held_atac_exposed"),
    "l40s_validation_reused": False,
    "science_shape_unchanged_vs_v10": True,
    "gpu_runs": len(grid),
    "cpu_control_runs": len(cpu_runs),
    "resource_firewall_snapshot_sha256": first["snapshot_sha256"],
    "resource_quiet_seconds_at_first_probe": round(quiet_seconds, 1),
    "resource_newest_protected_file": newest_path,
    "resource_probe_count_agreeing": 4,
    "cross_model_comparison_constraint": (
        "Cobolt v11 executes on gpu_rna_atac_torch_b6k, a runtime no other "
        "observed_multiome member shares. Its registered comparator is the "
        "150-run CPU control matrix re-fit under this campaign, and all 25 GPU "
        "fits run serially on one B6K device, so within-bundle hardware "
        "homogeneity holds by construction. A cross-model observed_multiome "
        "ranking may not treat this campaign as runtime-comparable with "
        "L40S-fit members without first establishing numeric interchangeability "
        "under config/b6k_numeric_tolerance_prespecification_20260825.json."
    ),
    "supersedes_cpu_control_matrix": (
        "executions/cobolt-v9-current-cpu-baseline-bundles/"
        "cobolt-v9-cpu-baselines--21088479 is superseded FOR v11 COMPARISON "
        "PURPOSES ONLY; it remains intact as its own frozen v9 record."
    ),
}
(stage / "receipt.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
freeze_tree(stage, {"artifact_class": "cobolt_v11_prefreeze_resource_stability"})
print(json.dumps(receipt, indent=2, sort_keys=True))
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Cobolt v11 frozen; audit at %s\n' "${TARGET}"
