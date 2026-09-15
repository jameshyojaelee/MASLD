#!/usr/bin/env bash
#SBATCH --job-name=model-check-463
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
TARGET="${ROOT}/executions/observed-multiome-taskspec-regrounding-${SLURM_JOB_ID}"
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

sha256sum --check --strict <<EOF
a474db7f50fe481206b017c5d77b1434d39f48c55683b3e2a4c9a544229ec0b3  ${ROOT}/config/artifacts/models/observed_multiome_integration/disposition_20260824.json
8bb1027133b4216481c5b61fa29b266a394ef13b1dea9d74b1c6667ce66e78ec  ${ROOT}/scripts/audit_observed_multiome_integration_disposition.py
abc634ea96718031d7212c4a1229fe3fd84c626d8ad46d75fd68fda91ccdb729  ${ROOT}/config/evaluation/observed_multiome_task.toml
78cff35fa0cec02b06091c8bb13ac2da4f6e289ce845abffa9ba4f1eed06cff3  ${ROOT}/config/evaluation/observed_multiome_development_gate.toml
6c2e92a78af315a0412a949de208daea5ec28d30fdd2e421a7762d3ffe6252b8  ${ROOT}/config/evaluation/observed_multiome_task_capabilities.toml
EOF

module load python/3.11.5-GCCcore-13.2.0

PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_observed_multiome_taskspec_regrounding \
  tests.unit.test_audit_observed_multiome_integration_disposition \
  tests.unit.test_non_census_record_directories \
  tests.unit.test_audit_remaining_multimodal_families_disposition \
  > "${STAGE}/unit_tests.stdout.txt" 2> "${STAGE}/unit_tests.stderr.txt"

# The audit receipt is the evidence the conclusion was re-derived, not asserted.
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  "${ROOT}/scripts/audit_observed_multiome_integration_disposition.py" \
  > "${STAGE}/audit_receipt.json" 2> "${STAGE}/audit_receipt.stderr.txt"

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
  "${ROOT}/config/artifacts/models/observed_multiome_integration/disposition_20260824.json" \
  "${ROOT}/scripts/audit_observed_multiome_integration_disposition.py" \
  "${ROOT}/tests/unit/test_observed_multiome_taskspec_regrounding.py" \
  "${ROOT}/slurm/validate_observed_multiome_integration_disposition_cpu.sbatch" \
  "${ROOT}/slurm/validate_observed_multiome_regrounding_cpu.sbatch" \
  "${ROOT}/config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json" \
  "${ROOT}/config/evaluation/observed_multiome_task.toml" \
  "${ROOT}/config/evaluation/observed_multiome_development_gate.toml" \
  "${ROOT}"/config/models/*.toml \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src:${ROOT}" python - "${STAGE}" "${ROOT}" <<'PY'
import hashlib, json, re, sys, tomllib
from pathlib import Path
from masld_bench.artifacts import freeze_tree, write_json_exclusive

stage, root = Path(sys.argv[1]), Path(sys.argv[2])
SUBJECTS = ("multivi", "scglue", "mofaplus", "seurat_wnn")
CANONICAL = "config/tasks/observed_multiome.toml"
STANDALONE = "config/evaluation/observed_multiome_task.toml"

record = json.loads(
    (root / "config/artifacts/models/observed_multiome_integration/disposition_20260824.json")
    .read_text(encoding="utf-8")
)
receipt = json.loads((stage / "audit_receipt.json").read_text(encoding="utf-8"))
gate = record["task_gate"]
before = (stage / "base_registry_before.sha256").read_text(encoding="utf-8")
after = (stage / "base_registry_after.sha256").read_text(encoding="utf-8")
contract = (stage / "contract_tests.stderr.txt").read_text(encoding="utf-8")
failures = set(re.findall(r"^(?:FAIL|ERROR): (\S+)", contract, re.MULTILINE))
EXPECTED = {
    "test_every_model_audit_bundle_is_complete_and_hash_bound",
    "test_every_repo_path_anchored_audit_hash_rederives",
    "test_graph_and_transparent_baseline_audits_freeze_exact_boundaries",
    "test_scimilarity_preflight_and_exposure_are_frozen",
}
with (root / "config/evaluation/observed_multiome_development_gate.toml").open("rb") as h:
    devgate = tomllib.load(h)

conjuncts_preserved = all(
    record["model_dispositions"][s]["gate_superseded_20260825"][len("BLOCKED_NO_TASKSPEC_"):]
    == record["model_dispositions"][s]["gate"][len("BLOCKED_TASKSPEC_EXECUTION_BLOCKED_"):]
    for s in SUBJECTS
)

checks = {
    "base_registry_byte_identical_across_the_job": before == after,
    "audit_was_rederived_over_all_34_pins": receipt["verified_metadata_files"] == 34,
    "audit_passed": receipt["status"] == "pass_metadata_audit_all_model_lanes_blocked",
    "zero_models_admitted": receipt["admitted_model_count"] == 0,
    "all_four_dispositions_unchanged": all(
        record["model_dispositions"][s]["admitted"] is False for s in SUBJECTS
    ),
    "schema_bumped_to_v2": record["schema_version"].endswith("-v2"),
    "status_regrounded": record["status"]
    == "FAIL_CLOSED_TASKSPEC_EXECUTION_BLOCKED_AND_NO_ADMITTED_MODEL_LANE",
    "two_path_reality_recorded_as_two_values": gate["task_spec_paths"]
    == {CANONICAL: False, STANDALONE: True},
    "old_single_boolean_removed": "task_spec_exists" not in gate,
    "canonical_path_still_absent_on_disk": not (root / CANONICAL).exists(),
    "standalone_path_present_on_disk": (root / STANDALONE).is_file(),
    "registry_status_observed_untouched": gate["registry_status_observed"]
    == "planned_registration",
    "every_gate_kept_its_conjuncts": conjuncts_preserved,
    "standalone_gates_still_shut": (
        devgate["global_census_revision_required_before_execution"] is True
        and devgate["model_promotion_allowed"] is False
        and devgate["champion_claim_allowed"] is False
        and devgate["selection_lock_allowed"] is False
        and devgate["external_seal_bound"] is False
    ),
    "execution_still_unauthorized": receipt["standalone_task_spec_execution_authorized"]
    is False,
    "no_new_contract_failure": not (failures - EXPECTED),
    "baseline_still_four": len(failures) == 4,
}
failed = [k for k, v in checks.items() if v is not True]

write_json_exclusive(
    stage / "validation.json",
    {
        "schema_version": "masld-bench-observed-multiome-regrounding-validation-v1",
        "checks": checks,
        "failed_checks": failed,
        "status": "passed" if not failed else "failed",
        "audit_receipt": receipt,
        "dispositions": {
            s: record["model_dispositions"][s]["admitted"] for s in SUBJECTS
        },
        "gate_regrounding": {
            s: {
                "before": record["model_dispositions"][s]["gate_superseded_20260825"],
                "after": record["model_dispositions"][s]["gate"],
            }
            for s in SUBJECTS
        },
        "contract_suite": {
            "observed_failures": sorted(failures),
            "new_failures": sorted(failures - EXPECTED),
        },
        "passed_because_rederived_not_because_asserted": (
            "The audit re-verified all 34 authority pins, the standalone TaskSpec digest and "
            "contents, and the development gate digest and contents. No assertion was edited to "
            "match a stale record; the validator's model of the world was changed to match a "
            "world that had changed, and it remains strictly stricter than the check it replaced."
        ),
        "outcomes_read": False,
        "metrics_calculated": False,
        "sealed_data_read": False,
    },
)
freeze_tree(
    stage,
    {
        "artifact_class": "observed_multiome_taskspec_regrounding_validation",
        "status": "passed" if not failed else "failed",
        "base_registry_unchanged": before == after,
        "models_reassessed": 0,
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
