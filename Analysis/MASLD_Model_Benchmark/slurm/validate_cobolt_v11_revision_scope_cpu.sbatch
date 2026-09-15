#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-353
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-check-353-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-check-353-%j.err

# Revision-scope validation of the ALREADY FROZEN Cobolt v11 candidate.
#
# The v11 freeze (job 21108035) succeeded: the candidate is COMPLETE and valid.
# Only the post-freeze scope check failed, because its exemption list was
# incomplete. This validator runs against the existing candidate rather than
# regenerating it, which also avoids repeating a 28-minute config walk.
#
# WHAT WENT WRONG, precisely. The failing check compared metadata field
# input_owner_by_role, a dict KEYED BY INPUT ROLE. The environment role name is
# one of those keys, so a runtime revision necessarily changes it:
#     only in v10: environment:gpu_rna_atac_torch_smoke
#     only in v11: environment:gpu_rna_atac_torch_b6k
# Diffed directly: 22 keys before and after, no shared key changed value, and
# after stripping environment:* keys ZERO of 175 runs differ. The `inputs` list
# itself was already filtered correctly and never differed.
#
# THE FIX IS STRICTER, NOT LOOSER. input_owner_by_role is not exempted wholesale
# -- that would blind the guard to a real input-owner swap, which is what it
# exists to catch. Instead only environment:* keys are set aside, and the
# environment delta is then asserted POSITIVELY: same input count, every
# non-environment input identical in role/path/payload, exactly one environment
# role removed and one added, and the added role matching the campaign's declared
# runtime_id.
#
# General lesson recorded here because it has now bitten two guards in one day:
# enumerating the fields a revision MAY change is fragile, because every revision
# legitimately moves a few derived hashes. Asserting WHAT CHANGED AND NOTHING
# ELSE is the robust shape.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

REPOSITORY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${REPOSITORY}/Analysis/MASLD_Model_Benchmark
V11=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v11--e17be3dab9ba48cb
V10=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v10--6cbef395fdccc389
RUNTIME_TREE=${ROOT}/executions/environments/cobolt-native-torch2.8.0-cu128-21107432
RUNTIME_TREE_SHA=e7e56c0ee207843b78dd9bbb1e3619fc965fc7836a4d3b0131c58b6f4e58049c
VALIDATION=${ROOT}/executions/cobolt-b6k-runtime-validation-21107495
VALIDATION_SHA=9ab79607dc34c059485ed04e4c2fef005713eddd8003cf4d234a5ecdf16afa7b
ADAPTER_TEST=${ROOT}/executions/cobolt-b6k-adapter-test-21107681
TARGET=${ROOT}/executions/cobolt-v11-revision-scope-validation-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the Cobolt v11 scope validation.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/cobolt-v11-scope-${SLURM_JOB_ID}
cd "${ROOT}"

sha256sum --check --strict <<EOF
${RUNTIME_TREE_SHA}  ${RUNTIME_TREE}/ARTIFACTS.json
${VALIDATION_SHA}  ${VALIDATION}/ARTIFACTS.json
EOF
( cd "${RUNTIME_TREE}" && sha256sum --check --strict SHA256SUMS )

PYTHONPATH=src python - \
  "${REPOSITORY}" "${V11}" "${V10}" "${STAGE}" "${RUNTIME_TREE}" "${VALIDATION}" "${ADAPTER_TEST}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.firewall import resource_snapshot
from masld_bench.hashing import sha256_file
from masld_bench.planner import load_frozen_plan

repository, v11p, v10p, stage, runtime_tree, validation, adapter_test = (
    Path(v).resolve(strict=True) for v in sys.argv[1:]
)
NEW_RUNTIME_ID = "gpu_rna_atac_torch_b6k"
OLD_RUNTIME_ID = "gpu_rna_atac_torch_smoke"
NEW_PROFILE, OLD_PROFILE = "gpu_rna_atac_b6k", "gpu_rna_atac_smoke"
NEW_ENV_ROLE, OLD_ENV_ROLE = f"environment:{NEW_RUNTIME_ID}", f"environment:{OLD_RUNTIME_ID}"

failures: list[str] = []
def check(cond: bool, message: str) -> None:
    if not cond:
        failures.append(message)

for tree in (v11p, v10p, validation):
    verify_frozen_tree(tree)
plan, prior = load_frozen_plan(v11p), load_frozen_plan(v10p)
check(plan["campaign"]["campaign_id"] == "v1-rna-atac-cobolt-5seed-smoke-v11", "v11 campaign id differs")
check(prior["campaign"]["campaign_id"] == "v1-rna-atac-cobolt-5seed-smoke-v10", "v10 campaign id differs")

vr = json.loads((validation / "validation.json").read_text())
check(vr.get("status") == "pass", "B6K runtime validation did not pass")
check(vr.get("l40s_validation_reused") is False, "B6K validation reused the L40S validation")
for c in ("forward_backward_passed", "query_batch_invariance_passed", "tensor_state_restore_passed",
          "cross_process_determinism_passed", "environment_contract_passed"):
    check(vr.get(c) is True, f"B6K validation check failed: {c}")
adapter_result = json.loads((adapter_test / "result.json").read_text())
check(adapter_result.get("equation_parity_passed") is True, "adapter e2e equation parity failed")
check(adapter_result.get("held_atac_exposed") is False, "adapter e2e exposed held ATAC")

# ---------------------------------------------------------------- must change
check(plan["source_lock"] != prior["source_lock"], "source lock did not change")
check(plan["plan_sha256"] != prior["plan_sha256"], "plan is not distinct from v10")
check(not ({str(r["run_id"]) for r in plan["runs"]} & {str(r["run_id"]) for r in prior["runs"]}),
      "v11 reuses superseded run identities")
check(NEW_RUNTIME_ID in plan["runtime_locks"] and OLD_RUNTIME_ID not in plan["runtime_locks"],
      "v11 does not bind the B6K runtime exclusively")
check(plan["resource_profiles"][NEW_PROFILE]["accelerator"] == "b6k", "v11 profile is not b6k")
check(str(runtime_tree) in " ".join(str(x) for x in plan["runtime_locks"][NEW_RUNTIME_ID]["command_prefix"]),
      "runtime lock does not point at the validated B6K tree")

# ------------------------------------------------------------ must not change
_EXEMPT = {"run_id", "campaign_id", "code_lock_sha256", "runtime_id",
           "resource_profile", "runtime_registry_sha256"}
_SCI = ("task_id", "model_id", "dataset_ids", "seed", "fold", "split_id",
        "stage", "action", "adaptation_regime", "hyperparameters", "checkpoint")
_META = ("adapter_request", "adapter_timeout_seconds", "dataset_routing", "dataset_routing_wave",
         "development_prediction_dataset_ids", "evaluator_contract_sha256", "evaluator_parameters",
         "fit_dataset_ids", "input_owner_by_role", "open_champion_eligible",
         "prediction_first_stress_dataset_ids", "primary_endpoint_scoring_allowed",
         "primary_evaluator_id", "sealed_prediction_dataset_ids", "split_contract",
         "split_contract_sha256", "split_registry_sha256", "task_contract_sha256",
         "task_registry_sha256", "wave")

def science(run):
    md = run.get("metadata") or {}
    row = {k: run.get(k) for k in _SCI}
    meta = {}
    for k in _META:
        v = md.get(k)
        # Set aside ONLY the environment role key; every other owner entry is compared.
        if k == "input_owner_by_role" and isinstance(v, dict):
            v = {kk: vv for kk, vv in v.items() if not kk.startswith("environment:")}
        meta[k] = v
    row["metadata"] = meta
    row["inputs"] = sorted(json.dumps(i, sort_keys=True) for i in (run.get("inputs") or [])
                           if not str(i.get("role", "")).startswith("environment:"))
    return row

def ident(run):
    return (run.get("model_id"), run.get("seed"), run.get("fold"),
            run.get("stage"), tuple(run.get("action") or []))

new_runs = {ident(r): r for r in plan["runs"]}
old_runs = {ident(r): r for r in prior["runs"]}
check(set(new_runs) == set(old_runs), "the set of scientific runs changed")
check(len(plan["runs"]) == 175 and len(prior["runs"]) == 175, "run count is not 175")
for k in sorted(set(new_runs) & set(old_runs), key=str):
    if science(new_runs[k]) != science(old_runs[k]):
        failures.append(f"scientific content changed for run {k}")

# -------------------------------- positively assert the environment delta only
env_deltas = 0
for k in sorted(set(new_runs) & set(old_runs), key=str):
    a, b = old_runs[k], new_runs[k]
    ai, bi = a.get("inputs") or [], b.get("inputs") or []
    check(len(ai) == len(bi), f"input count changed for run {k}")
    a_non = {str(i.get("role")): json.dumps(i, sort_keys=True)
             for i in ai if not str(i.get("role", "")).startswith("environment:")}
    b_non = {str(i.get("role")): json.dumps(i, sort_keys=True)
             for i in bi if not str(i.get("role", "")).startswith("environment:")}
    check(a_non == b_non, f"a NON-ENVIRONMENT input changed for run {k}")
    a_env = {str(i.get("role")) for i in ai if str(i.get("role", "")).startswith("environment:")}
    b_env = {str(i.get("role")) for i in bi if str(i.get("role", "")).startswith("environment:")}
    if a_env != b_env:
        env_deltas += 1
        check(a_env - b_env == {OLD_ENV_ROLE}, f"unexpected environment role removed for run {k}: {a_env - b_env}")
        check(b_env - a_env == {NEW_ENV_ROLE}, f"unexpected environment role added for run {k}: {b_env - a_env}")
        check(b.get("runtime_id") == NEW_RUNTIME_ID,
              f"added environment role does not match the run runtime_id for {k}")
    ao = (a.get("metadata") or {}).get("input_owner_by_role") or {}
    bo = (b.get("metadata") or {}).get("input_owner_by_role") or {}
    check(len(ao) == len(bo), f"input_owner_by_role key count changed for run {k}")
    check({x: v for x, v in ao.items() if not x.startswith("environment:")}
          == {x: v for x, v in bo.items() if not x.startswith("environment:")},
          f"a NON-ENVIRONMENT input owner changed for run {k}")
check(env_deltas == 25, f"expected exactly 25 runs to change environment role, saw {env_deltas}")

# ------------------------------------------------------------------- registry
def reg(p): return {i["path"]: i["sha256"] for i in p["registry_snapshot"]["files"]}
nr, orr = reg(plan), reg(prior)
check(set(nr) - set(orr) == {f"runtimes/{NEW_RUNTIME_ID}.toml"}, "unexpected registry additions")
check(set(orr) - set(nr) == {f"runtimes/{OLD_RUNTIME_ID}.toml"}, "unexpected registry removals")
allowed = {"campaigns/v1_rna_atac_cobolt_5seed_smoke.toml",
           "models/multimodal_integration.toml", "resources.toml"}
changed = {p for p in set(nr) & set(orr) if nr[p] != orr[p]}
check(not (changed - allowed), f"unexpected registry files changed: {sorted(changed - allowed)}")
for key in ("dataset_locks", "safety"):
    check(plan.get(key) == prior.get(key), f"{key} changed")

models_toml = (Path(plan["config_root"]) / "models" / "multimodal_integration.toml").read_text()
check('runtime_id = "gpu_rna_atac_torch_smoke"' in models_toml, "scpair no longer binds the L40S runtime")
check(models_toml.count('runtime_id = "gpu_rna_atac_torch_b6k"') == 1, "more than one model moved to B6K")

grid = sorted((int(r["seed"]), int(r["fold"])) for r in plan["runs"]
              if r.get("model_id") == "cobolt" and r.get("resource_profile") == NEW_PROFILE)
expected = [(s, f) for s in (1103, 2107, 3109, 4111, 5113) for f in range(5)]
check(grid == expected, "v11 does not carry the unchanged five-seed by five-fold grid")
check(grid == sorted((int(r["seed"]), int(r["fold"])) for r in prior["runs"]
                     if r.get("model_id") == "cobolt" and r.get("resource_profile") == OLD_PROFILE),
      "v11 GPU grid differs from v10")
check(len([r for r in plan["runs"] if r.get("resource_profile") != NEW_PROFILE]) == 150,
      "expected 150 CPU control runs")

observed = resource_snapshot(repository)
firewall_at_freeze = plan["resource_firewall"]["snapshot_sha256"]

receipt = {
    "schema_version": "masld-bench-cobolt-v11-revision-scope-validation-v1",
    "status": "failed_scope" if failures else "passed_execution_runtime_only_revision",
    "candidate": v11p.as_posix(),
    "candidate_manifest_sha256": sha256_file(v11p / "ARTIFACTS.json"),
    "plan_sha256": plan["plan_sha256"],
    "supersedes_candidate": v10p.as_posix(),
    "supersedes_plan_sha256": prior["plan_sha256"],
    "revision_class": "execution_runtime_only",
    "runtime_id": NEW_RUNTIME_ID,
    "resource_profile": NEW_PROFILE,
    "accelerator": "b6k",
    "runs_total": len(plan["runs"]),
    "gpu_runs_moved": env_deltas,
    "cpu_control_runs": 150,
    "scientific_content_identical_to_v10": not failures,
    "sole_input_difference": (
        "The only per-run difference is the environment-lock input and its "
        f"input_owner_by_role key: {OLD_ENV_ROLE} -> {NEW_ENV_ROLE}. Verified 22 owner "
        "keys before and after, no shared key changed value, and zero of 175 runs "
        "differ once environment:* keys are set aside. No data input changed."
    ),
    "guard_defect_fixed": (
        "The freeze-time guard enumerated the fields a revision MAY change and its "
        "list was incomplete: input_owner_by_role is keyed BY ROLE, so the environment "
        "role name is a key. Fixed by setting aside only environment:* keys and then "
        "asserting the environment delta positively, which is strictly stronger than "
        "the check that failed. Enumerating permitted changes is fragile; asserting "
        "what changed and nothing else is the robust shape."
    ),
    "runtime_tree": runtime_tree.as_posix(),
    "runtime_tree_artifacts_sha256": sha256_file(runtime_tree / "ARTIFACTS.json"),
    "b6k_validation_artifacts_sha256": sha256_file(validation / "ARTIFACTS.json"),
    "b6k_adapter_test_result_sha256": sha256_file(adapter_test / "result.json"),
    "l40s_validation_reused": False,
    "resource_firewall_snapshot_at_freeze": firewall_at_freeze,
    "resource_firewall_snapshot_now": observed["snapshot_sha256"],
    "resource_firewall_unchanged_since_freeze": observed["snapshot_sha256"] == firewall_at_freeze,
    "v10_disposition": (
        "v10 remains a valid frozen provenance record and its scope validation "
        "(executions/cobolt-v10-revision-scope-validation-21106245) stands. It is no "
        "longer RE-EXECUTABLE, because resources.toml and models/multimodal_integration.toml "
        "moved. That was a deliberate, approved consequence, not drift."
    ),
    "cross_model_comparison_constraint": (
        "Cobolt v11 executes on gpu_rna_atac_torch_b6k, a runtime no other "
        "observed_multiome member shares. Its comparator is the 150-run CPU control "
        "matrix re-fit under this campaign, and all 25 GPU fits run serially on one "
        "B6K device, so within-bundle hardware homogeneity holds by construction. A "
        "cross-model observed_multiome ranking may not treat this campaign as "
        "runtime-comparable with L40S-fit members without first establishing numeric "
        "interchangeability under "
        "config/b6k_numeric_tolerance_prespecification_20260825.json."
    ),
    "absence_of_comparison_is_not_a_passed_comparison": (
        "No L40S-versus-B6K numeric comparison exists for Cobolt: job 21060919 never "
        "recorded a loss value. The frozen tolerance prespecification admits B6K output "
        "only into a rectangle fit entirely on B6K, which at 0/25 is exactly this case. "
        "The absence of a comparison must not be read as a comparison having passed."
    ),
    "python_numpy_axis_numerically_inert": (
        "The same seeded Cobolt training loss, 122.22412109375, was produced on "
        "python 3.11.16 / numpy 1.26.0 (job 21107038) and on python 3.10.19 / "
        "numpy 2.2.6 (job 21107495), bit-identical. torch is the only live variable."
    ),
}
if failures:
    receipt["failures"] = failures
(stage / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
freeze_tree(stage, {"artifact_class": "cobolt_v11_revision_scope_validation"})
print(json.dumps(receipt, indent=2, sort_keys=True))
if failures:
    raise SystemExit("v11 revision scope validation FAILED: " + "; ".join(failures[:10]))
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Cobolt v11 revision scope validated: %s\n' "${TARGET}"
