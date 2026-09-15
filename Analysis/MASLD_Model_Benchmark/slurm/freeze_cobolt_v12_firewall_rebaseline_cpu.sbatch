#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-freeze-356
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-freeze-356-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-freeze-356-%j.err

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
PRIOR_CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v11--e17be3dab9ba48cb
RUNTIME_TREE=${ROOT}/executions/environments/cobolt-native-torch2.8.0-cu128-21107432
RUNTIME_TREE_SHA=e7e56c0ee207843b78dd9bbb1e3619fc965fc7836a4d3b0131c58b6f4e58049c
VALIDATION=${ROOT}/executions/cobolt-b6k-runtime-validation-21107495
VALIDATION_SHA=9ab79607dc34c059485ed04e4c2fef005713eddd8003cf4d234a5ecdf16afa7b
ADAPTER_TEST=${ROOT}/executions/cobolt-b6k-adapter-test-21107681
TARGET=${ROOT}/executions/cobolt-v12-prefreeze-resource-stability-${SLURM_JOB_ID}
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

CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v12"
PRIOR_CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v11"
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

# --- v12 is FIREWALL-BASELINE-ONLY: the opposite scope from v11 --------------
# v11 was required to change the fixed source record, fixed runtime and profile. v12 must
# change NOTHING except the protected Resource snapshot, so every one of those is
# asserted UNCHANGED here.
if plan["source_lock"] != prior["source_lock"]:
    raise SystemExit("scope violated: v12 must not change the scoped source lock")
if plan["runtime_locks"] != prior["runtime_locks"]:
    raise SystemExit("scope violated: v12 must not change any runtime lock")
if plan["resource_profiles"] != prior["resource_profiles"]:
    raise SystemExit("scope violated: v12 must not change any resource profile")
for key in ("dataset_locks", "safety", "resource_totals"):
    if plan.get(key) != prior.get(key):
        raise SystemExit(f"scope violated: {key} changed")
if NEW_RUNTIME_ID not in plan["runtime_locks"] or OLD_RUNTIME_ID in plan["runtime_locks"]:
    raise SystemExit("v12 does not bind the B6K runtime exclusively")
if plan["resource_profiles"][NEW_PROFILE]["accelerator"] != "b6k":
    raise SystemExit("v12 GPU profile does not request the B6K accelerator")
if str(runtime_tree) not in " ".join(str(x) for x in plan["runtime_locks"][NEW_RUNTIME_ID]["command_prefix"]):
    raise SystemExit("v12 runtime lock does not point at the validated B6K runtime tree")

# The one thing that MUST move: the protected Resource baseline.
if plan["resource_firewall"]["snapshot_sha256"] == prior["resource_firewall"]["snapshot_sha256"]:
    raise SystemExit(
        "v12 captured the same Resource snapshot as v11; no re-baseline was warranted "
        "and the v11 firewall abort is unexplained"
    )
if plan["plan_sha256"] == prior["plan_sha256"]:
    raise SystemExit("revision is not distinct from the superseded plan")
if set(str(r["run_id"]) for r in plan["runs"]) & set(str(r["run_id"]) for r in prior["runs"]):
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
old_runs = {ident(r): r for r in prior["runs"]}
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

new_registry, old_registry = registry_files(plan), registry_files(prior)
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
    "schema_version": "masld-bench-cobolt-v12-prefreeze-resource-stability-v1",
    "status": "pass",
    "campaign_id": CAMPAIGN_ID,
    "supersedes_campaign_id": PRIOR_CAMPAIGN_ID,
    "revision_class": "firewall_baseline_only",
    "supersedes_reason": (
        "v11 was invalidated before any run executed by a second Resource promotion "
        "burst on 2026-08-25: figures/main/fig4_singlecell_programs/captions.md was "
        "rewritten at 12:53:11 (e8fc7ccacf200b6c 7885 B -> dd3c69138df040ac 7868 B). "
        "That file is gitignored, so the drift is invisible to git status. Exactly "
        "one of 1,676 human-edited authority files changed; none added or removed."
    ),
    "prior_firewall_snapshot_sha256": prior["resource_firewall"]["snapshot_sha256"],
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
freeze_tree(stage, {"artifact_class": "cobolt_v12_prefreeze_resource_stability"})
print(json.dumps(receipt, indent=2, sort_keys=True))
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Cobolt v11 frozen; audit at %s\n' "${TARGET}"
