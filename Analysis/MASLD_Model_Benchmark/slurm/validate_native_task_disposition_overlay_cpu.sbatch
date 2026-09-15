#!/usr/bin/env bash
#SBATCH --job-name=model-check-461
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
TARGET="${ROOT}/executions/native-task-disposition-overlay-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"; FAILED="${TARGET}.failed"
preserve() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve EXIT
[[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]] && { echo refusing >&2; exit 1; }
install -d -m 0750 "${STAGE}"
cd "${ROOT}"

# ---------------------------------------------------------------------------
# 1. The base registry is the surface this overlay must not touch.  Hash it
#    before anything runs, and hash it again at the end.
# ---------------------------------------------------------------------------
BASE_BEFORE="${STAGE}/base_registry_before.sha256"
BASE_AFTER="${STAGE}/base_registry_after.sha256"
sha256sum "${ROOT}"/config/models/*.toml > "${BASE_BEFORE}"

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
# 2. The overlay itself: 22 record-only dispositions plus the test and this
#    job script.
# ---------------------------------------------------------------------------
sha256sum --check --strict <<EOF
46c9d44af22e1ff054bf33d442a8c9178f7d08dc36ff20e1b4f1a4680c8403e7  ${ROOT}/config/artifacts/models/alphagenome/native_task_disposition_20260825.json
b7d6ea23237d87fa8b698a69004dcea6c0c6913b7d538253e365688c5a00696f  ${ROOT}/config/artifacts/models/borzoi/native_task_disposition_20260825.json
25db07f2c65519d722c172459445a85cd1f19b3afb9ef9c82133f8f862d80323  ${ROOT}/config/artifacts/models/caduceus/native_task_disposition_20260825.json
169ce9e876ff17f64faf54c6e04af9f79414c9a82382f998e4d6c7d1edf8d376  ${ROOT}/config/artifacts/models/cellplm/native_task_disposition_20260825.json
a104d9f38ae74dbdb69a328cca471e253406fe74e384966698f566f1f4ae0147  ${ROOT}/config/artifacts/models/dnabert2/native_task_disposition_20260825.json
fab2dd781351ec935c627f8e11a19a30346c9682de5b448097b3097d0f6bf445  ${ROOT}/config/artifacts/models/enformer/native_task_disposition_20260825.json
a8b91ec9ccc3d08a95d7ad618f4155162cb5c0175a75213ac7b1103b83313ace  ${ROOT}/config/artifacts/models/epibert/native_task_disposition_20260825.json
714d42c80f945752df4da7b399bc6a7502ec8c081164a11c7fcd9b6118eeb4e9  ${ROOT}/config/artifacts/models/evo/native_task_disposition_20260825.json
299575830c9127b33b94e669f18db7d063da9bbedf0bae608b43ec9711b43a8b  ${ROOT}/config/artifacts/models/geneformer/native_task_disposition_20260825.json
b0f944059fa4c4a862c4a721d2ec85adcc5380f64565bee334aeb105c0ed9ac9  ${ROOT}/config/artifacts/models/hyenadna/native_task_disposition_20260825.json
73020610a36b2a0c5564a95014aef581912b8a8023f22df36f2762f71c1f4283  ${ROOT}/config/artifacts/models/langcell/native_task_disposition_20260825.json
341923c8ec73d122c239c77505e9187c2b934a7d66d5f41a39ddbc35bfbedcae  ${ROOT}/config/artifacts/models/nucleotide_transformer/native_task_disposition_20260825.json
7525f32822948ed724797c32c40272c79b71b115dd6e65567287946365d7c4a2  ${ROOT}/config/artifacts/models/regformer/native_task_disposition_20260825.json
b263886aa8a2809c2ce54ea88703bd5ae2b4de24603b9cc9d039ac473417d0a3  ${ROOT}/config/artifacts/models/scbert/native_task_disposition_20260825.json
233f50ea7603a83745e6764b2bf8508c7ad6c29b50b1f7f48f78540fc2032acd  ${ROOT}/config/artifacts/models/sccello/native_task_disposition_20260825.json
19146c55928c3a44e8ae20ca04d11ac42642678d14664071e09d8bfcc41a3dad  ${ROOT}/config/artifacts/models/scgpt/native_task_disposition_20260825.json
e06e5580fef691647d3c6bb7490eda7898233c31f213b0841cc87b75b20af428  ${ROOT}/config/artifacts/models/scimilarity/native_task_disposition_20260825.json
ce98bb5131830aeb7308df1defa4769d7263d6ad935313ca5d5579017bf32ba4  ${ROOT}/config/artifacts/models/scooby_epicardioids/native_task_disposition_20260825.json
0f37532bfaea4f87b48fdc33ac6f6ed123d3f804ecf55fc4ea6bcf9548e86451  ${ROOT}/config/artifacts/models/scooby_neurips/native_task_disposition_20260825.json
8c5ce1295017cda060cda2453e0455c3178c0f729c7f59c6ff673e8ee1cd1ee6  ${ROOT}/config/artifacts/models/scooby_onek1k/native_task_disposition_20260825.json
d359e662b44e160fe082e9ab83d7f3a33d678cae26eb1dc1441b424e9314de68  ${ROOT}/config/artifacts/models/sei/native_task_disposition_20260825.json
a9525e62c2567418b4eec17d7ee510feb32bf240da3cdc32517ed06c8e0645ee  ${ROOT}/config/artifacts/models/uce/native_task_disposition_20260825.json
003d17b02f756216bc8965d4084eef0b9e973c5c97c0080546cac3eca4677bca  ${ROOT}/tests/unit/test_native_task_disposition_overlay.py
EOF

module load python/3.11.5-GCCcore-13.2.0

# ---------------------------------------------------------------------------
# 3. Overlay unit tests.
# ---------------------------------------------------------------------------
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_native_task_disposition_overlay \
  > "${STAGE}/unit_tests.stdout.txt" 2> "${STAGE}/unit_tests.stderr.txt"

# ---------------------------------------------------------------------------
# 4. Full contract suite.  Nine tests in tests/contract/test_registry_files.py
#    were already failing before this overlay existed, for reasons this record
#    does not touch: stale config/models/*.toml source pins inside other audit
#    bundles, two audit-bundle directories with no registry model, and
#    scgpt/uce/scimilarity preflight assertions that disagree with their own
#    bundles.  The check here is therefore a regression check: the failing set
#    must be exactly the pre-existing set, no more.
# ---------------------------------------------------------------------------
set +e
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest discover -s tests/contract -t . -v \
  > "${STAGE}/contract_tests.stdout.txt" 2> "${STAGE}/contract_tests.stderr.txt"
CONTRACT_RC=$?
set -e
echo "${CONTRACT_RC}" > "${STAGE}/contract_tests.returncode.txt"

python3 -m pip freeze --all > "${STAGE}/pip_freeze.txt" 2>/dev/null || true

sha256sum "${ROOT}"/config/models/*.toml > "${BASE_AFTER}"

sha256sum \
  "${ROOT}"/config/artifacts/models/*/native_task_disposition_20260825.json \
  "${ROOT}/tests/unit/test_native_task_disposition_overlay.py" \
  "${ROOT}/slurm/validate_native_task_disposition_overlay_cpu.sbatch" \
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
NAME = "native_task_disposition_20260825.json"
OVERLAY = "model_census_native_task_disposition_20260825"
MODELS = root / "config" / "artifacts" / "models"

# Failures that already existed before this overlay was authored, confirmed on
# an isolated mirror of the tree with the 22 records removed.
PREEXISTING = {
    "test_all_registry_files_parse_and_have_expected_schema",
    "test_every_campaign_binds_the_frozen_model_census",
    "test_every_model_audit_bundle_is_complete_and_hash_bound",
    "test_every_repo_path_anchored_audit_hash_rederives",
    "test_graph_and_transparent_baseline_audits_freeze_exact_boundaries",
    "test_model_census_maps_exactly_to_audit_bundles",
    "test_scgpt_upstream_preflight_and_exposure_are_frozen",
    "test_scimilarity_preflight_and_exposure_are_frozen",
    "test_uce_upstream_preflight_and_exposure_are_frozen",
}

records = {
    path.parent.name: json.loads(path.read_text(encoding="utf-8"))
    for path in sorted(MODELS.glob(f"*/{NAME}"))
}

# --- base registry byte-identity across the whole job ----------------------
before = (stage / "base_registry_before.sha256").read_text(encoding="utf-8")
after = (stage / "base_registry_after.sha256").read_text(encoding="utf-8")

# --- every bound evidence digest re-derives from the live file -------------
evidence_checked = 0
evidence_bad = []
for family, record in records.items():
    for item in record["bound_evidence"]:
        target = root / item["path"]
        actual = hashlib.sha256(target.read_bytes()).hexdigest()
        evidence_checked += 1
        if actual != item["sha256"]:
            evidence_bad.append(f"{family}: {item['path']}")

# --- the base registry rows still say what each disposition was read from --
rows = {}
for path in sorted((root / "config" / "models").glob("*.toml")):
    with path.open("rb") as handle:
        for model in tomllib.load(handle).get("models", []):
            rows[model["model_id"]] = model

drift = []
covered = set()
for family, record in records.items():
    for model_id, facts in record["registry_facts_at_audit"].items():
        covered.add(model_id)
        row = rows[model_id]
        for key in (
            "license_status",
            "exposure_status",
            "checkpoint_sha256",
            "status",
            "admission_blocking",
        ):
            if row[key] != facts[key]:
                drift.append(f"{model_id}:{key}")
        if list(row["supported_tasks"]) != facts["supported_tasks"]:
            drift.append(f"{model_id}:supported_tasks")

by_task = {}
for record in records.values():
    by_task.setdefault(record["task_id"], set()).update(record["model_ids"])

# --- requirements-suite regression check ---------------------------------------
contract = (stage / "contract_tests.stderr.txt").read_text(encoding="utf-8")
observed_failures = set(
    re.findall(r"^(?:FAIL|ERROR): (\S+)", contract, re.MULTILINE)
)
new_failures = sorted(observed_failures - PREEXISTING)
resolved = sorted(PREEXISTING - observed_failures)
ran = re.search(r"^Ran (\d+) tests", contract, re.MULTILINE)

checks = {
    "twenty_two_records_present": len(records) == 22,
    "every_record_is_record_only": all(
        r["record_only"] and not r["consumed_by_runtime"] for r in records.values()
    ),
    "every_record_declares_base_registry_unchanged": all(
        r["base_registry_unchanged"] for r in records.values()
    ),
    "every_record_carries_the_overlay_id": all(
        r["overlay_id"] == OVERLAY for r in records.values()
    ),
    "no_record_read_an_outcome": all(
        not r["outcomes_read"] and not r["sealed_data_read"]
        for r in records.values()
    ),
    "no_record_calculated_a_metric": all(
        not r["metrics_calculated"] for r in records.values()
    ),
    "no_record_opened_a_checkpoint": all(
        not r["checkpoint_opened"]
        and not r["unrestricted_deserializer_used"]
        and not r["license_gate_bypassed"]
        for r in records.values()
    ),
    "no_record_authorizes_scoring": all(
        r["binding_consequences"]["sealed_scoring_authorized"] is False
        for r in records.values()
    ),
    "no_record_authorizes_a_champion": all(
        r["binding_consequences"]["sealed_champion_eligible"] is False
        and r["binding_consequences"]["open_champion_eligible"] is False
        for r in records.values()
    ),
    "base_registry_byte_identical_across_the_job": before == after,
    "base_registry_digests_match_the_records": all(
        hashlib.sha256((root / item["registry_file"]).read_bytes()).hexdigest()
        == item["registry_file_sha256"]
        for r in records.values()
        for item in r["base_registry_files"]
    ),
    "every_bound_evidence_digest_rederives": not evidence_bad,
    "registry_facts_have_not_drifted": not drift,
    "cell_state_mapping_covers_fourteen_model_ids": len(
        by_task.get("cell_state_mapping", ())
    ) == 14,
    "sequence_native_regulatory_covers_ten_model_ids": len(
        by_task.get("sequence_native_regulatory", ())
    ) == 10,
    "observed_multiome_covers_three_model_ids": len(
        by_task.get("observed_multiome", ())
    ) == 3,
    "twenty_seven_model_ids_covered_in_total": len(covered) == 27,
    "overlay_unit_tests_passed": True,
    "no_new_contract_failure": not new_failures,
    "contract_failure_count_is_unchanged": len(observed_failures) == len(PREEXISTING),
}

failed = [key for key, value in checks.items() if value is not True]

write_json_exclusive(
    stage / "validation.json",
    {
        "schema_version": "masld-bench-native-task-disposition-overlay-validation-v1",
        "overlay_id": OVERLAY,
        "checks": checks,
        "failed_checks": failed,
        "status": "passed" if not failed else "failed",
        "bound_evidence_digests_checked": evidence_checked,
        "bound_evidence_mismatches": evidence_bad,
        "registry_fact_drift": drift,
        "model_ids_by_task": {
            task: sorted(ids) for task, ids in sorted(by_task.items())
        },
        "contract_suite": {
            "tests_run": int(ran.group(1)) if ran else None,
            "returncode": int(
                (stage / "contract_tests.returncode.txt").read_text().strip()
            ),
            "preexisting_failures": sorted(PREEXISTING),
            "observed_failures": sorted(observed_failures),
            "new_failures_introduced_by_this_overlay": new_failures,
            "previously_failing_now_passing": resolved,
            "note": (
                "The nine pre-existing failures were confirmed on an isolated "
                "mirror of the tree with all 22 overlay records removed. They "
                "are stale config/models/*.toml source pins inside other audit "
                "bundles, two audit-bundle directories with no registry model, "
                "and scgpt/uce/scimilarity preflight assertions that disagree "
                "with their own bundles. This overlay changes none of them."
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
        "artifact_class": "native_task_disposition_overlay_validation",
        "overlay_id": OVERLAY,
        "status": "passed" if not failed else "failed",
        "records": len(records),
        "model_ids_covered": len(covered),
        "base_registry_unchanged": before == after,
        "new_contract_failures": len(new_failures),
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
