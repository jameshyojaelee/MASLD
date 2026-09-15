#!/usr/bin/env bash
#SBATCH --job-name=model-check-464
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/scimilarity-corrections-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"; FAILED="${TARGET}.failed"
preserve() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve EXIT
[[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]] && { echo refusing >&2; exit 1; }
install -d -m 0750 "${STAGE}"
cd "${ROOT}"

sha256sum "${ROOT}"/config/models/*.toml > "${STAGE}/base_registry_before.sha256"
sha256sum --check --strict <<EOF
6c0d8d6aff2db6fa99000ff5cc593940f097b6dde9745d08b8cd93aed9b4d14d  ${ROOT}/config/models/cell_foundation.toml
5fabd7dafffa3c2233a752d8a5016c11670b115bae2a0e4e5f3fc6a8123984ab  ${ROOT}/config/models/chromatin.toml
bb6e11d6b38b2a424f719a867d92e88c5d9e275223cf06385e82c89f6621e386  ${ROOT}/config/models/dna_language.toml
926641e4df3746e5497b91f8518565b259ebc7fa41837f36f5e021d5ee83ae10  ${ROOT}/config/models/evidence_graph.toml
cd7e67aa068c2a041fbf842564fc20d82be008020bb5226d70bc36c937e4a74d  ${ROOT}/config/models/mandatory_baselines.toml
9569f89cf5d452e0335f899c3a3e43a483661fa2259d471d2208a5b1c5c107cb  ${ROOT}/config/models/multimodal_integration.toml
30afc3d28d6ecae4f3924eee0a149e0a55859514ee9c5dbb54c79c847109e10d  ${ROOT}/config/models/perturbation.toml
8c6684a23ef40e6573c788fdfc9d697668f5f5fc978390118dc8f12f9c79ae1d  ${ROOT}/config/models/regulatory_local.toml
bccce82a8a949eb81b73246848dcf7e80b164e464992fdb9b9fcb10c72181c1b  ${ROOT}/config/models/regulatory_sequence.toml
1bbe758b19e13989c795ec1d1c444523ee4a50e89ea99823194d5458fad4c364  ${ROOT}/config/models/restricted_unclear.toml
e611592a542fd4d205aeaa79d89dcefa5aa2b783cc4b197bae97c4ff8ae8c426  ${ROOT}/config/models/spatial_deferred.toml
EOF

# The three scimilarity bundle files are UNCHANGED by this campaign. Only the
# test's accepted sets and two assertions moved. Pinned here to prove it.
sha256sum --check --strict <<EOF
73d6a2083d763e9afd1eee2276a1df08119ebf9e40ebaf4115e95d32b59f44e8  ${ROOT}/config/artifacts/models/scimilarity/checkpoints.json
142dcf3e57361e81e43f1fd8fd715d776a915ba43b4bd0dd047b76aa6b247eca  ${ROOT}/config/artifacts/models/scimilarity/exposure_audit.json
1688aa2b8a535bde4f3b3c2de11fc9fb0403117657e79913d61a99c719a24709  ${ROOT}/config/artifacts/models/scimilarity/development_crosswalk.json
9fb381fdbdc2a3758f40d7b877cf66ab35c0555d373594bd23a75c4df9ddcdd7  ${ROOT}/tests/contract/test_registry_files.py
faea633576c079c1d8db50173dc0dad40544a9b067eff5b5ce4f54b0ea3f1fcb  ${ROOT}/tests/unit/test_scimilarity_admission_and_schema_correction.py
1c58a57d742a7a254c63f5ff863dff6cb232a0c03a96afd40dbdb22608d3ced4  ${ROOT}/scripts/check_contract_suite_baseline.py
fd711d8d3f3aa792e2aa7c076bf43c41126acb4f9a27edcc202762be5474521f  ${ROOT}/config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json
EOF

module load python/3.11.5-GCCcore-13.2.0

PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_scimilarity_admission_and_schema_correction \
  tests.unit.test_observed_multiome_taskspec_regrounding \
  tests.unit.test_non_census_record_directories \
  tests.unit.test_native_task_disposition_overlay \
  > "${STAGE}/unit_tests.stdout.txt" 2> "${STAGE}/unit_tests.stderr.txt"

set +e
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest discover -s tests/contract -t . -v \
  > "${STAGE}/contract_tests.stdout.txt" 2> "${STAGE}/contract_tests.stderr.txt"
echo "$?" > "${STAGE}/contract_tests.returncode.txt"
set -e

PYTHONPATH="${ROOT}/src:${ROOT}" python "${ROOT}/scripts/check_contract_suite_baseline.py" \
  --unittest-output "${STAGE}/contract_tests.stderr.txt" \
  --require-module tests.contract.test_registry_files \
  --minimum-tests 125 \
  > "${STAGE}/baseline_gate.stdout.txt" 2> "${STAGE}/baseline_gate.stderr.txt"

python3 -m pip freeze --all > "${STAGE}/pip_freeze.txt" 2>/dev/null || true
sha256sum "${ROOT}"/config/models/*.toml > "${STAGE}/base_registry_after.sha256"
sha256sum \
  "${ROOT}/tests/contract/test_registry_files.py" \
  "${ROOT}/tests/unit/test_scimilarity_admission_and_schema_correction.py" \
  "${ROOT}/scripts/check_contract_suite_baseline.py" \
  "${ROOT}/config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json" \
  "${ROOT}/slurm/validate_scimilarity_corrections_cpu.sbatch" \
  "${ROOT}"/config/artifacts/models/scimilarity/*.json \
  "${ROOT}"/config/models/*.toml \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src:${ROOT}" python - "${STAGE}" "${ROOT}" <<'PY'
import hashlib, json, re, sys, tomllib
from pathlib import Path
from masld_bench.artifacts import freeze_tree, write_json_exclusive

stage, root = Path(sys.argv[1]), Path(sys.argv[2])
MODELS = root / "config" / "artifacts" / "models"
SC = MODELS / "scimilarity"

models = {}
for path in sorted((root / "config" / "models").glob("*.toml")):
    with path.open("rb") as handle:
        for m in tomllib.load(handle).get("models", []):
            models[m["model_id"]] = m
sci = models["scimilarity_v1_1"]
violations = [
    mid for mid, m in models.items() if m["admission_blocking"] != bool(m["blockers"])
]

def shape(node, prefix=""):
    out = set()
    if isinstance(node, dict):
        for k, v in node.items():
            out.add(prefix + k); out |= shape(v, f"{prefix}{k}.")
    elif isinstance(node, list) and node and isinstance(node[0], dict):
        out |= shape(node[0], f"{prefix}[].")
    return out

superset = {}
for fname, stem in (
    ("checkpoints.json", "masld-bench-upstream-checkpoint-preflight-v"),
    ("development_crosswalk.json", "masld-bench-development-crosswalk-v"),
    ("exposure_audit.json", "masld-bench-checkpoint-exposure-audit-v"),
):
    peers = []
    for b in sorted(p for p in MODELS.iterdir() if p.is_dir()):
        m = b / fname
        if not m.is_file():
            continue
        body = json.loads(m.read_text(encoding="utf-8"))
        if body.get("schema_version") == f"{stem}1":
            peers.append(shape(body))
    v2 = shape(json.loads((SC / fname).read_text(encoding="utf-8")))
    superset[fname] = {
        "v1_peers": len(peers),
        "universal_paths_missing_from_v2": sorted(set.intersection(*peers) - v2),
    }

before = (stage / "base_registry_before.sha256").read_text(encoding="utf-8")
after = (stage / "base_registry_after.sha256").read_text(encoding="utf-8")
contract = (stage / "contract_tests.stderr.txt").read_text(encoding="utf-8")
failures = set(re.findall(r"^(?:FAIL|ERROR): (\S+)", contract, re.MULTILINE))
ran = re.search(r"^Ran (\d+) tests", contract, re.MULTILINE)
EXPECTED = {
    "test_every_repo_path_anchored_audit_hash_rederives",
    "test_graph_and_transparent_baseline_audits_freeze_exact_boundaries",
}

checks = {
    "base_registry_byte_identical_across_the_job": before == after,
    "admission_blocking_invariant_holds_for_all_136": not violations and len(models) == 136,
    "scimilarity_still_admission_blocked": sci["admission_blocking"] is True,
    "scimilarity_blocker_names_the_contract_not_the_exposure_class": (
        "ModelExecutionContract" in sci["blockers"][0]
        and "Keep this model admission-blocked" in sci["blockers"][0]
    ),
    "scimilarity_still_the_only_encoder_seen_model": sorted(
        m for m, v in models.items() if v["exposure_status"] == "encoder_seen"
    ) == ["scimilarity_v1_1"],
    "no_registry_value_was_changed": before == after,
    "every_v2_file_is_a_superset_of_the_v1_shape": all(
        not v["universal_paths_missing_from_v2"] for v in superset.values()
    ),
    "both_scimilarity_tests_now_pass": not (
        {
            "test_every_model_audit_bundle_is_complete_and_hash_bound",
            "test_scimilarity_preflight_and_exposure_are_frozen",
        }
        & failures
    ),
    "no_new_contract_failure": not (failures - EXPECTED),
    "baseline_is_two": len(failures) == 2,
    "only_pin_cascade_remains": failures == EXPECTED,
}
failed = [k for k, v in checks.items() if v is not True]

write_json_exclusive(
    stage / "validation.json",
    {
        "schema_version": "masld-bench-scimilarity-corrections-validation-v1",
        "checks": checks,
        "failed_checks": failed,
        "status": "passed" if not failed else "failed",
        "item_1_schema_evolution": {
            "verdict": "sanctioned_evolution_test_was_behind",
            "superset_check": superset,
            "fix": "widened the accepted sets; bundle untouched",
        },
        "item_2_admission_blocking": {
            "verdict": "registry_correct_test_was_wrong",
            "invariant_violations_across_the_census": violations,
            "models_checked": len(models),
            "scimilarity_admission_blocking": sci["admission_blocking"],
            "scimilarity_blocker_count": len(sci["blockers"]),
            "fix": "corrected the two assertions to the invariant; registry untouched",
            "two_gates": (
                "The admission blocker is the missing ModelExecutionContract. encoder_seen drives "
                "sealed-champion ineligibility, a different gate. Clearing one would not clear the "
                "other."
            ),
        },
        "contract_suite": {
            "tests_run": int(ran.group(1)) if ran else None,
            "observed_failures": sorted(failures),
            "resolved_by_this_campaign": [
                "test_every_model_audit_bundle_is_complete_and_hash_bound",
                "test_scimilarity_preflight_and_exposure_are_frozen",
            ],
            "new_failures": sorted(failures - EXPECTED),
        },
        "outcomes_read": False,
        "metrics_calculated": False,
        "sealed_data_read": False,
    },
)
freeze_tree(
    stage,
    {
        "artifact_class": "scimilarity_corrections_validation",
        "status": "passed" if not failed else "failed",
        "base_registry_unchanged": before == after,
        "bundles_changed": 0,
        "contract_failures_now": len(failures),
        "outcomes_read": False,
        "metrics_calculated": False,
        "sealed_data_read": False,
    },
)
print(json.dumps({"status": "passed" if not failed else "failed", "failed": failed}))
if failed:
    raise SystemExit("validation failed: " + ", ".join(failed))
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
