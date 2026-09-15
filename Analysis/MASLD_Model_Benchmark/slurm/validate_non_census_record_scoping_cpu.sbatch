#!/usr/bin/env bash
#SBATCH --job-name=model-check-462
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
TARGET="${ROOT}/executions/non-census-record-scoping-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"; FAILED="${TARGET}.failed"
preserve() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve EXIT
[[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]] && { echo refusing >&2; exit 1; }
install -d -m 0750 "${STAGE}"
cd "${ROOT}"

# ---------------------------------------------------------------------------
# 1. The frozen base registry is the surface this campaign must not touch.
#    Hash it first, hash it again at the end, and assert equality in-job.
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# 2. Everything this campaign authored or refreshed, including the full pinner
#    cascade for both records so nothing is left stranded.
# ---------------------------------------------------------------------------
sha256sum --check --strict <<EOF
4cebd852fa7a4b64ee139d8a264f637c29f98198d9696c5ffeb028185b295f2a  ${ROOT}/config/evaluation/model_audit_bindings.toml
9fb381fdbdc2a3758f40d7b877cf66ab35c0555d373594bd23a75c4df9ddcdd7  ${ROOT}/tests/contract/test_registry_files.py
3ef6162a4b764714c02f2272200b9f3e27019a8e5cae3f9d834319b2b0f988fe  ${ROOT}/tests/unit/test_non_census_record_directories.py
fd711d8d3f3aa792e2aa7c076bf43c41126acb4f9a27edcc202762be5474521f  ${ROOT}/config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json
1c58a57d742a7a254c63f5ff863dff6cb232a0c03a96afd40dbdb22608d3ced4  ${ROOT}/scripts/check_contract_suite_baseline.py
a474db7f50fe481206b017c5d77b1434d39f48c55683b3e2a4c9a544229ec0b3  ${ROOT}/config/artifacts/models/observed_multiome_integration/disposition_20260824.json
edf10ff7709c83d3db36389e2bdeea301ecc11637a53b038fc7799cc932c2cad  ${ROOT}/slurm/validate_observed_multiome_integration_disposition_cpu.sbatch
de07b243f55b38a616fee8c333540a50baa2ec284f2d08944b6f5a259171e591  ${ROOT}/config/artifacts/models/remaining_multimodal_families/disposition_20260824.json
e2a93d052b40fe3d452d9e2e702b627c56829b3524ed1cbceb579fce1fc2a12b  ${ROOT}/config/artifacts/models/remaining_multimodal_families/validation_21088037.json
6ab00be39fde06d411e2ac8836ef298555927c0da4746331d99e313da713f471  ${ROOT}/slurm/validate_remaining_multimodal_families_disposition_cpu.sbatch
960c81e7132bceeff473966a8554bc56522f15e561c65ba7a1e4c0cd500d3ddf  ${ROOT}/slurm/validate_native_task_disposition_overlay_cpu.sbatch
EOF

module load python/3.11.5-GCCcore-13.2.0

# ---------------------------------------------------------------------------
# 3. Unit tests: the anti-laundering guard, the two refreshed record
#    validators, and the earlier native-task overlay whose pins were refreshed.
#    remaining_multimodal_families must pass. observed_multiome_integration is
#    expected to FAIL on a substantive re-audit demand its stale pin was
#    hiding, so it is recorded rather than conditional; see the debt register.
# ---------------------------------------------------------------------------
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_non_census_record_directories \
  tests.unit.test_native_task_disposition_overlay \
  tests.unit.test_audit_remaining_multimodal_families_disposition \
  > "${STAGE}/unit_tests.stdout.txt" 2> "${STAGE}/unit_tests.stderr.txt"

set +e
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_audit_observed_multiome_integration_disposition \
  > "${STAGE}/observed_multiome_record_audit.stdout.txt" \
  2> "${STAGE}/observed_multiome_record_audit.stderr.txt"
echo "$?" > "${STAGE}/observed_multiome_record_audit.returncode.txt"
set -e

# ---------------------------------------------------------------------------
# 4. Full contract suite through the project's own baseline check.
# ---------------------------------------------------------------------------
set +e
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest discover -s tests/contract -t . -v \
  > "${STAGE}/contract_tests.stdout.txt" 2> "${STAGE}/contract_tests.stderr.txt"
echo "$?" > "${STAGE}/contract_tests.returncode.txt"
set -e

PYTHONPATH="${ROOT}/src:${ROOT}" python "${ROOT}/scripts/check_contract_suite_baseline.py" \
  --unittest-output "${STAGE}/contract_tests.stderr.txt" \
  --require-module tests.contract.test_registry_files \
  --minimum-tests 105 \
  > "${STAGE}/baseline_gate.stdout.txt" 2> "${STAGE}/baseline_gate.stderr.txt"

python3 -m pip freeze --all > "${STAGE}/pip_freeze.txt" 2>/dev/null || true
sha256sum "${ROOT}"/config/models/*.toml > "${STAGE}/base_registry_after.sha256"

sha256sum \
  "${ROOT}/config/evaluation/model_audit_bindings.toml" \
  "${ROOT}/tests/contract/test_registry_files.py" \
  "${ROOT}/tests/unit/test_non_census_record_directories.py" \
  "${ROOT}/config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json" \
  "${ROOT}/scripts/check_contract_suite_baseline.py" \
  "${ROOT}/slurm/validate_non_census_record_scoping_cpu.sbatch" \
  "${ROOT}"/config/artifacts/models/observed_multiome_integration/*.json \
  "${ROOT}"/config/artifacts/models/remaining_multimodal_families/*.json \
  "${ROOT}"/config/models/*.toml \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src:${ROOT}" python - "${STAGE}" "${ROOT}" <<'PY'
import hashlib
import json
import re
import sys
import tomllib
from pathlib import Path

from masld_bench.artifacts import freeze_tree, write_json_exclusive

stage, root = Path(sys.argv[1]), Path(sys.argv[2])
AUDIT = root / "config" / "artifacts" / "models"
REQUIRED = {"checkpoints.json", "development_crosswalk.json", "exposure_audit.json"}

authority = tomllib.loads(
    (root / "config/evaluation/model_audit_bindings.toml").read_text(encoding="utf-8")
)
declared = set(authority["non_census_record_directories"])
rules = authority["non_census_record_rules"]

models = {}
for path in sorted((root / "config" / "models").glob("*.toml")):
    with path.open("rb") as handle:
        for model in tomllib.load(handle).get("models", []):
            models[model["model_id"]] = model
grouped = {
    m: b["bundle_id"] for b in authority["group_bindings"] for m in b["model_ids"]
}
resolved = {m: grouped.get(m, m) for m in models}
bundle_ids = set(resolved.values())
observed = {p.name for p in AUDIT.iterdir() if p.is_dir()}

# The proof by construction: registering the two records instead of declaring
# them re-resolves their subjects away from their own directories.
def orphans_if_registered(extra):
    g = dict(grouped)
    for bid, subjects in extra:
        for m in subjects:
            g[m] = bid
    ids = {m: g.get(m, m) for m in models}
    return len(observed - set(ids.values())), len(set(ids.values()))

def named(directory):
    text = "".join(p.read_text(encoding="utf-8") for p in sorted(directory.glob("*.json")))
    return sorted(
        m for m in models
        if re.search(r"(?<![A-Za-z0-9_])" + re.escape(m) + r"(?![A-Za-z0-9_])", text)
    )

subjects = {name: named(AUDIT / name) for name in sorted(declared)}
registered_orphans, registered_bundles = orphans_if_registered(
    [(n, s) for n, s in subjects.items()]
)

before = (stage / "base_registry_before.sha256").read_text(encoding="utf-8")
after = (stage / "base_registry_after.sha256").read_text(encoding="utf-8")
contract = (stage / "contract_tests.stderr.txt").read_text(encoding="utf-8")
observed_failures = set(re.findall(r"^(?:FAIL|ERROR): (\S+)", contract, re.MULTILINE))
ran = re.search(r"^Ran (\d+) tests", contract, re.MULTILINE)

RESOLVED_HERE = "test_model_census_maps_exactly_to_audit_bundles"
EXPECTED_STILL_FAILING = {
    "test_every_model_audit_bundle_is_complete_and_hash_bound",
    "test_every_repo_path_anchored_audit_hash_rederives",
    "test_graph_and_transparent_baseline_audits_freeze_exact_boundaries",
    "test_scimilarity_preflight_and_exposure_are_frozen",
}

checks = {
    "base_registry_byte_identical_across_the_job": before == after,
    "two_directories_declared": sorted(declared)
    == ["observed_multiome_integration", "remaining_multimodal_families"],
    "declaration_does_not_shadow_a_census_bundle": not (declared & bundle_ids),
    "audit_root_is_exactly_census_plus_declared": observed == bundle_ids | declared,
    "declared_dirs_hold_no_census_member_file": all(
        not ({p.name for p in (AUDIT / n).iterdir() if p.is_file()} & REQUIRED)
        for n in declared
    ),
    "declared_dirs_hold_a_disposition_record": all(
        sorted((AUDIT / n).glob(rules["disposition_record_glob"])) for n in declared
    ),
    "every_named_subject_owns_its_own_bundle": all(
        resolved[m] == m and (AUDIT / m).is_dir()
        for names in subjects.values()
        for m in names
    ),
    "aggregate_rule_still_has_zero_exceptions": all(
        not (AUDIT / m).is_dir()
        for b in authority["group_bindings"]
        for m in b["model_ids"]
    ),
    "registering_them_would_be_strictly_worse": registered_orphans > len(declared),
    "census_mapping_test_now_passes": RESOLVED_HERE not in observed_failures,
    "no_new_contract_failure": not (observed_failures - EXPECTED_STILL_FAILING),
    "baseline_shrank_from_five_to_four": len(observed_failures) == 4,
    "remaining_multimodal_families_validator_passes": True,
    "base_registry_model_count_unchanged": len(models)
    == authority["expected_registry_model_count"],
    "expected_artifact_bundle_count_unchanged": len(bundle_ids)
    == authority["expected_artifact_bundle_count"],
}
failed = [k for k, v in checks.items() if v is not True]

write_json_exclusive(
    stage / "validation.json",
    {
        "schema_version": "masld-bench-non-census-record-scoping-validation-v1",
        "checks": checks,
        "failed_checks": failed,
        "status": "passed" if not failed else "failed",
        "declared_non_census_directories": sorted(declared),
        "declared_subjects": subjects,
        "census_bundle_count": len(bundle_ids),
        "audit_root_directory_count": len(observed),
        "proof_by_construction": {
            "as_declared": {"bundles": len(bundle_ids), "orphans": len(declared)},
            "if_registered_as_group_bindings": {
                "bundles": registered_bundles,
                "orphans": registered_orphans,
            },
            "conclusion": (
                "Registering the two records re-resolves their subjects away from "
                "their own directories and multiplies the orphan set, so they "
                "cannot be census members."
            ),
        },
        "contract_suite": {
            "tests_run": int(ran.group(1)) if ran else None,
            "observed_failures": sorted(observed_failures),
            "resolved_by_this_campaign": [RESOLVED_HERE],
            "still_failing_by_design": sorted(EXPECTED_STILL_FAILING),
            "new_failures": sorted(observed_failures - EXPECTED_STILL_FAILING),
        },
        "observed_multiome_record_audit": {
            "returncode": int(
                (stage / "observed_multiome_record_audit.returncode.txt").read_text().strip()
            ),
            "expected_to_fail": True,
            "reason": (
                "Refreshing its stale multimodal_integration.toml pin unmasked a "
                "substantive staleness: the record asserts task_spec_exists=false but "
                "config/evaluation/observed_multiome_task.toml was created 2026-08-25. "
                "Its validator correctly demands a new audit, which only the "
                "observed-multiome lane can perform. Not closed here."
            ),
        },
        "outcomes_read": False,
        "metrics_calculated": False,
        "sealed_data_read": False,
    },
)

freeze_tree(
    stage,
    {
        "artifact_class": "non_census_record_scoping_validation",
        "status": "passed" if not failed else "failed",
        "base_registry_unchanged": before == after,
        "contract_failures_now": len(observed_failures),
        "resolved_by_this_campaign": 1,
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
