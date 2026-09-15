#!/usr/bin/env python3
"""Reconcile remaining cell specialists and one donor-safe 50k execution route."""

from __future__ import annotations

import argparse
import hashlib
import json
import tomllib
from pathlib import Path
from typing import Any

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


class CellSpecialistAuditError(RuntimeError):
    """Raised when a cell-specialist authority or disposition differs."""


CENSUS_RELATIVE = (
    "config/artifacts/model-authorities/remaining_cell_specialists/census.json"
)
CENSUS_SHA256 = "cfea4be9c745dc2a061a004fceb5178d35c467cafec373e8ca02dfc534f077ea"

AUTHORITY_SHA256 = {
    "config/artifacts/models/scgpt/checkpoints.json": "41333904c75f9c513ad74d113a4dd7a1b97300e3cb71e3ac442c56a8cdd0943f",
    "config/artifacts/models/scgpt/development_crosswalk.json": "d6a9e3f182423778027afad116aa20806a001c2a480f74238e8ccbb3c408699a",
    "config/artifacts/models/scgpt/exposure_audit.json": "b0078da8341965fedfa2508d840c70c1ea5e13ea4866fbf5ca5242582629f9b7",
    "config/artifacts/models/scimilarity/checkpoints.json": "73d6a2083d763e9afd1eee2276a1df08119ebf9e40ebaf4115e95d32b59f44e8",
    "config/artifacts/models/scimilarity/development_crosswalk.json": "1688aa2b8a535bde4f3b3c2de11fc9fb0403117657e79913d61a99c719a24709",
    "config/artifacts/models/scimilarity/exposure_audit.json": "142dcf3e57361e81e43f1fd8fd715d776a915ba43b4bd0dd047b76aa6b247eca",
    "config/artifacts/models/cellplm/checkpoints.json": "eecb7f6a82541edc9524bfa3fa5c1c06731acb88cdfe6025a08204e4defb4fca",
    "config/artifacts/models/cellplm/development_crosswalk.json": "b5156f78e91d18ff443f2aad27ab4d1877708d549529eaacae335ef6e67de865",
    "config/artifacts/models/cellplm/exposure_audit.json": "8e1bae0d482d8c4b51977ffb7abb78819b3271f3d4c1edf52cf0cde47a9e7626",
    "config/artifacts/models/scprint/checkpoints.json": "ebcf61926e8af397e32ab8ad4e26f1618c53595e4309efcbf485ef3ba9ac3b1a",
    "config/artifacts/models/scprint/development_crosswalk.json": "479ba64fe4e45b1fb154724d151e836d85b5ad05163ba8362163b8d8fd4de624",
    "config/artifacts/models/scprint/exposure_audit.json": "bbf35d2a0c0f05397beb5d2924165679582992ea049291733118138fcebc8bf9",
    "config/artifacts/models/scprint/runtime_disposition.json": "66df641faf6505d4992a4633aa94e015dd1187e495972c16cfed8bc0b07fce4e",
    "config/artifacts/models/scprint2/checkpoints.json": "14a979f9f0afef6822ae3d9575c0007db4e110dc3a0cffd0452b2a05e02a1b96",
    "config/artifacts/models/scprint2/development_crosswalk.json": "2ac7927f82804a0140c39368c382f0b47d2c7cc6361a73d43414c7ee32e7853b",
    "config/artifacts/models/scprint2/exposure_audit.json": "4222198d978b737d6653009da4ac25b346a21bd40a507a2cd4e275f39ab3f364",
    "config/artifacts/models/regformer/checkpoints.json": "9886bde7b8804ea25ff7d4800a4d418047cd0de10e3926a8c2609de6596cfdf7",
    "config/artifacts/models/regformer/development_crosswalk.json": "4d046904751aed4c1ff89ba47ef585adb0bae824365ee00448e2429bfc047b3f",
    "config/artifacts/models/regformer/exposure_audit.json": "9a8d9bd2743e1ba0c72c98f109a878f5096737396625ac90b5d972cfd28fb763",
    "config/artifacts/models/scbert/checkpoints.json": "a35dbfdd8bfa5b890eb17209e25ca3fb4c82fe05631758d536958df3d7836deb",
    "config/artifacts/models/scbert/development_crosswalk.json": "8f31c7ee5667f929f9c8908932ed89ab113a9acd9c87b9a5a0925b5874c011ca",
    "config/artifacts/models/scbert/exposure_audit.json": "3d15f9ffaf8cb1141e45c90b2d3bae45b49b993bde842787433aa6fc0b829caf",
    "config/artifacts/models/cell_bulk_classical_baselines/checkpoints.json": "886468e16d92e8df0249b79a978d00239f278f6d59fbb60f54ad637ba1f5f938",
    "config/artifacts/models/cell_bulk_classical_baselines/development_crosswalk.json": "e5b6a5c5369dd363a386f75bd5a1b3a4b8c8f09e73eb865d1a09c39363c0b465",
    "config/artifacts/models/cell_bulk_classical_baselines/exposure_audit.json": "0d95f120b4ca1ba9db92e365e04504b8c642f2980d47ee8a2b884b0f569abffc",
}

CONFIG_SHA256 = {
    "config/models/cell_foundation.toml": "6c0d8d6aff2db6fa99000ff5cc593940f097b6dde9745d08b8cd93aed9b4d14d",
    "config/models/mandatory_baselines.toml": "cd7e67aa068c2a041fbf842564fc20d82be008020bb5226d70bc36c937e4a74d",
    "config/tasks/cell_state_mapping.toml": "2a07a6e5558b43d9825b7f3d5581ab418432320e31acff4067a9aab6dd9f678a",
}

EVIDENCE_SHA256 = {
    "executions/scgpt_whole_human-weights-only-inspection-21066250/weights_only_inspection.json": "3967e3f9211da72d26efb302ed94916caf5a237b61e47d5ca35e6a5458b19028",
    "executions/scgpt_continual-weights-only-inspection-21066249/weights_only_inspection.json": "555852cb46a4918b248eb0e4eabab86ff278c6608b33d84d359e25c33332697d",
    "executions/regformer-weights-only-inspection-21066098/weights_only_inspection.json": "28217cb286650ca4dd1e767870ef8a5e6496cc63e44410cfde87170ad2b459df",
    "executions/scprint_v1_5_medium-weights-only-inspection-21066255/weights_only_inspection.json": "a518630cfdb6e3291645ccf40e709802954e11920d14f5b8056c3d6afe02f59d",
    "scripts/extract_scimilarity_v1_1_frozen_screen_50000_resumable.py": "e664cef093641202dd43989bd284e3c126ec644262913a39f6bc5d67ad0b15bd",
    "slurm/run_scimilarity_v1_1_frozen_screen_50000_resumable.sbatch": "e35a7f5b5c6d7bbce6ffea53bda3fcd237db87658ee5d43bbd01fbb9938fc878",
}

ARTIFACT_SHA256 = {
    "executions/scimilarity-v1-1-compat-runtime-21079070": "8f66e2e478742665072a23ed4b115346830fcee516094c59f20164c94b3a5a4a",
    "executions/scimilarity-v1-1-frozen-screen-fixture-21079292": "b6d50557107b5e6a878480741ec6a4c781663218876943ae205578b1a4ecd1e5",
    "executions/scimilarity-v1-1-frozen-screen-preflight-21079643": "8dcf5c8d91f59ebfa5cbb9d55be7d3299c0eeaae6522b9be7c3a893306eeb6e9",
    "executions/model-training-504-queue-validation-v1": "df1b5cdc490da90379524f305bbc6357580f32ec063570ceaffea785025c581e",
    "executions/scprint-v1-runtime-preflight-21072885": "f0ed73826f18a07ebc42dc1e80ff4ff41475b5546d5b3da8dad8bcfb015eec23",
    "executions/environments/scprint2-native-torch2.8.0-v2-21073597": "4474a0d069d09a9a805e7d6e33936af95eeab9fb72d6e3c131e7e7733350842d",
    "executions/scprint2-native-runtime-preflight-v2-21074027": "9ed6d6db96bca280b2b2a99f8c3256be329978ad145e234aaac4815ef33295fa",
    "executions/scprint2-no-neighbor-fixture-21075645": "edf71c212340ebfdcc266ef7c0ecb73449193b36a17147e53d3801804082c9f7",
    "executions/model-work-302-21075963": "3c57bcf5f9c821e395966809033b4319971432102714c0051f2c8d449e654183",
    "executions/model-evaluation-302-21076583": "109077b2d8c90870946c406933ab2018a08db79e1965bcfad06198deb8b8245e",
}

DEEP_VERIFY = {
    "executions/scimilarity-v1-1-frozen-screen-preflight-21079643",
    "executions/model-training-504-queue-validation-v1",
    "executions/model-work-302-21075963",
    "executions/model-evaluation-302-21076583",
}

LEGACY_CHECKSUM_VERIFY = {
    "executions/scprint-v1-runtime-preflight-21072885": {
        "ARTIFACTS.json": "f0ed73826f18a07ebc42dc1e80ff4ff41475b5546d5b3da8dad8bcfb015eec23",
        "conda_explicit.txt": "e86f64c5d3670e841f4c586e5dbbe37603962fd9fb4c6b2c58fbcf9bc00b1002",
        "pip_freeze.txt": "5b2e1e6a63d1edf0b0d385194333d2fbd3adee793b9a9010d0cfffb02c51a9c0",
        "runtime_preflight.json": "ba264c18cda94353734cce59fa854353a53451064b3df26f641d49396be6beee",
    },
    "executions/scprint2-native-runtime-preflight-v2-21074027": {
        "ARTIFACTS.json": "9ed6d6db96bca280b2b2a99f8c3256be329978ad145e234aaac4815ef33295fa",
        "runtime_preflight.json": "a57677dac67329cd554b78f4e6dae51cfe4747bd32a940af5c1a0925606956d0",
    },
    "executions/scprint2-no-neighbor-fixture-21075645": {
        "ARTIFACTS.json": "edf71c212340ebfdcc266ef7c0ecb73449193b36a17147e53d3801804082c9f7",
        "fixture/ARTIFACTS.json": "5dc71a57ff3b18ff9d527639a253dc26a0eba7cca9cdf4023668448ebf793e06",
        "fixture/aligned_counts.h5ad": "951b2615c437c86991bf039b885c453c51d9af3b3d39a111cf33a914d8489635",
        "fixture/fixture_receipt.json": "ec842afe75a696205a4758a1d3dd6c56e12a93831a7ca7132aeaa19e8e721e5c",
    },
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CellSpecialistAuditError(f"JSON authority is unreadable: {path}: {error}") from error
    if not isinstance(value, dict):
        raise CellSpecialistAuditError(f"expected JSON object: {path}")
    return value


def load_toml(path: Path) -> dict[str, Any]:
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise CellSpecialistAuditError(f"TOML authority is unreadable: {path}: {error}") from error
    if not isinstance(value, dict):
        raise CellSpecialistAuditError(f"expected TOML table: {path}")
    return value


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CellSpecialistAuditError(message)


def require_hash(root: Path, relative: str, expected: str) -> Path:
    path = root / relative
    if not path.is_file() or digest(path) != expected:
        raise CellSpecialistAuditError(f"authority identity differs: {relative}")
    return path


def require_artifact(root: Path, relative: str, expected: str) -> dict[str, Any]:
    path = root / relative
    manifest_path = path / "ARTIFACTS.json"
    if not manifest_path.is_file() or digest(manifest_path) != expected:
        raise CellSpecialistAuditError(f"artifact manifest identity differs: {relative}")
    if relative in DEEP_VERIFY:
        try:
            return dict(verify_frozen_tree(path))
        except ArtifactError as error:
            raise CellSpecialistAuditError(f"frozen artifact differs: {relative}: {error}") from error
    if relative in LEGACY_CHECKSUM_VERIFY:
        expected_files = LEGACY_CHECKSUM_VERIFY[relative]
        observed_files = {
            child.relative_to(path).as_posix()
            for child in path.rglob("*")
            if child.is_file() and child.relative_to(path).as_posix() != "SHA256SUMS"
        }
        require(observed_files == set(expected_files),
                f"legacy artifact file census differs: {relative}")
        require(not any(child.is_symlink() for child in path.rglob("*")),
                f"legacy artifact contains a symlink: {relative}")
        for name, expected_file_sha256 in expected_files.items():
            require(digest(path / name) == expected_file_sha256,
                    f"legacy artifact payload differs: {relative}/{name}")
        ledger = (path / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
        expected_ledger = [
            f"{checksum}  {name}" for name, checksum in expected_files.items()
        ]
        require(ledger == expected_ledger,
                f"legacy artifact checksum ledger differs: {relative}")
    return load_json(manifest_path)


def model_rows(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    rows = config.get("models")
    if not isinstance(rows, list):
        raise CellSpecialistAuditError("model registry rows differ")
    indexed = {
        str(row.get("model_id")): row for row in rows if isinstance(row, dict)
    }
    if len(indexed) != len(rows):
        raise CellSpecialistAuditError("model registry contains duplicate or malformed IDs")
    return indexed


def validate_census(census: dict[str, Any]) -> tuple[set[str], dict[str, dict[str, Any]]]:
    expected = {
        "scgpt_whole_human",
        "scgpt_continual",
        "scimilarity_v1_1",
        "cellplm_85m",
        "cellplm_85m_vae_20231027",
        "scprint_v1_5_medium",
        "scprint2_small_v2",
        "scprint2_medium",
        "regformer",
        "scbert",
        "celltypist_baseline",
    }
    rows = census.get("models")
    require(isinstance(rows, list), "census model rows differ")
    indexed = {
        str(row.get("model_id")): row for row in rows if isinstance(row, dict)
    }
    require(set(indexed) == expected, "census model membership differs")
    require(len(indexed) == len(rows), "duplicate census model IDs")
    ready = [
        model_id
        for model_id, row in indexed.items()
        if row.get("direct_preregistered_wrapper_ready") is True
    ]
    require(ready == ["scimilarity_v1_1"], "direct-wrapper readiness differs")
    require(all(row.get("schedulable") is False for row in indexed.values()),
            "reconciliation silently authorized scheduling")
    require(all(row.get("sealed_champion_eligible") is False for row in indexed.values()),
            "reconciliation silently authorized a sealed champion")
    next_screen = census.get("next_executable_screen", {})
    require(next_screen.get("model_id") == "scimilarity_v1_1"
            and next_screen.get("adapter_fits_head") is False
            and next_screen.get("adapter_scores_model") is False
            and next_screen.get("submission_authorized") is False,
            "next 50k screen boundary differs")
    common = census.get("common_lane_contract", {})
    require(common.get("sealed_labels_read") is False
            and common.get("query_label_refinement_prohibited") is True
            and common.get("cross_donor_neighbor_or_context_coupling_prohibited") is True,
            "common-lane leakage boundary differs")
    require(census.get("scope", {}).get("sealed_sources_accessed") == [],
            "sealed source access was recorded")
    return expected, indexed


def audit(root: Path) -> dict[str, Any]:
    census = load_json(require_hash(root, CENSUS_RELATIVE, CENSUS_SHA256))
    expected, census_models = validate_census(census)
    authorities = {
        relative: load_json(require_hash(root, relative, expected_sha))
        for relative, expected_sha in AUTHORITY_SHA256.items()
    }
    configs = {
        relative: load_toml(require_hash(root, relative, expected_sha))
        for relative, expected_sha in CONFIG_SHA256.items()
    }
    evidence = {
        relative: load_json(require_hash(root, relative, expected_sha))
        if relative.endswith(".json") else {"sha256": expected_sha}
        for relative, expected_sha in EVIDENCE_SHA256.items()
    }
    artifacts = {
        relative: require_artifact(root, relative, expected_sha)
        for relative, expected_sha in ARTIFACT_SHA256.items()
    }

    foundation = model_rows(configs["config/models/cell_foundation.toml"])
    mandatory = model_rows(configs["config/models/mandatory_baselines.toml"])
    for model_id in expected - {"celltypist_baseline"}:
        row = foundation.get(model_id, {})
        require(row.get("supported_tasks") == ["cell_state_mapping"],
                f"{model_id} task scope differs")
        require(row.get("admission_blocking") is True,
                f"{model_id} was admitted without a new execution census")
    celltypist = mandatory.get("celltypist_baseline", {})
    require(celltypist.get("implementation_type") == "train_locally"
            and celltypist.get("checkpoint_sha256") == "NOT_APPLICABLE"
            and celltypist.get("admission_blocking") is True,
            "CellTypist registry boundary differs")
    task = configs["config/tasks/cell_state_mapping.toml"]
    require(task.get("unit_of_inference") == "donor", "cell task inferential unit differs")
    require("celltypist_baseline" in task.get("baseline_model_ids", []),
            "CellTypist is no longer a mandatory task baseline")

    scgpt = authorities["config/artifacts/models/scgpt/checkpoints.json"]
    scgpt_exposure = authorities["config/artifacts/models/scgpt/exposure_audit.json"]
    whole_inspection = evidence[
        "executions/scgpt_whole_human-weights-only-inspection-21066250/weights_only_inspection.json"
    ]
    continual_inspection = evidence[
        "executions/scgpt_continual-weights-only-inspection-21066249/weights_only_inspection.json"
    ]
    require(scgpt["bundles"]["scgpt_whole_human"]["artifacts"][0]["sha256"]
            == "UNRESOLVED", "scGPT registry checkpoint changed without a frozen bundle")
    require(whole_inspection.get("sha256")
            == census_models["scgpt_whole_human"]["observed_checkpoint_sha256"]
            and whole_inspection.get("tensor_count") == 163,
            "scGPT whole-human safe inventory differs")
    require(continual_inspection.get("sha256")
            == census_models["scgpt_continual"]["observed_checkpoint_sha256"]
            and continual_inspection.get("tensor_count") == 173,
            "scGPT continual safe inventory differs")
    require(scgpt.get("weight_license") == "UNDECLARED"
            and scgpt_exposure["checkpoint_findings"]["scgpt_continual"]["exposure_state"]
            == "target_label_unexposed", "scGPT terms or exposure differs")

    scimilarity = authorities["config/artifacts/models/scimilarity/checkpoints.json"]
    scimilarity_exposure = authorities[
        "config/artifacts/models/scimilarity/exposure_audit.json"
    ]
    require(scimilarity["staged_model_contract"]["encoder_sha256"]
            == census_models["scimilarity_v1_1"]["observed_checkpoint_sha256"],
            "SCimilarity encoder identity differs")
    require(scimilarity.get("code_license") == "Apache-2.0"
            and scimilarity.get("weight_license") == "CC-BY-SA-4.0",
            "SCimilarity terms differ")
    require(scimilarity_exposure["checkpoint_findings"]["scimilarity_v1_1"]["exposure_state"]
            == "encoder_seen", "SCimilarity development exposure differs")
    scim_runtime = artifacts["executions/scimilarity-v1-1-compat-runtime-21079070"]["metadata"]
    scim_fixture = artifacts["executions/scimilarity-v1-1-frozen-screen-fixture-21079292"]["metadata"]
    scim_preflight = artifacts["executions/scimilarity-v1-1-frozen-screen-preflight-21079643"]["metadata"]
    scim_queue = artifacts["executions/model-training-504-queue-validation-v1"]["metadata"]
    require(scim_runtime.get("status") == "passed_cpu_gate"
            and scim_runtime.get("strict_weights_only_restore") is True
            and scim_runtime.get("released_reference_index_used") is False,
            "SCimilarity runtime evidence differs")
    require(scim_fixture.get("development_rows") == 50000
            and scim_fixture.get("donors") == 102
            and scim_fixture.get("studies") == 7
            and scim_fixture.get("aggregate_exposure_status") == "encoder_seen",
            "SCimilarity 50k fixture differs")
    require(scim_preflight.get("status") == "passed"
            and scim_preflight.get("full_50000_extraction_submitted") is False
            and scim_preflight.get("sealed_champion_eligible") is False,
            "SCimilarity preflight boundary differs")
    require(scim_queue.get("status") == "passed"
            and scim_queue.get("gpu_submitted") is False
            and scim_queue.get("live_queue_mutated") is False
            and scim_queue.get("wrapper_sha256")
            == EVIDENCE_SHA256["slurm/run_scimilarity_v1_1_frozen_screen_50000_resumable.sbatch"],
            "SCimilarity queue contract differs")

    cellplm = authorities["config/artifacts/models/cellplm/checkpoints.json"]
    require(cellplm["terms_audit"]["result"] == "blocked_checkpoint_terms_undeclared"
            and cellplm["documented_checkpoint"]["sha256"] == "UNRESOLVED"
            and cellplm["registered_checkpoint"]["sha256"] == "UNRESOLVED",
            "CellPLM fail-closed release split differs")

    scprint_disposition = authorities[
        "config/artifacts/models/scprint/runtime_disposition.json"
    ]
    require(scprint_disposition.get("disposition") == "terminal_runtime_incompatibility"
            and scprint_disposition.get("strict_state_restoration_passed") is False
            and scprint_disposition.get("outcomes_or_labels_read") is False,
            "scPRINT v1.5 terminal disposition differs")

    scprint2 = authorities["config/artifacts/models/scprint2/checkpoints.json"]
    scprint2_exposure = authorities["config/artifacts/models/scprint2/exposure_audit.json"]
    require(scprint2["runtime_contract"]["environment_status"].startswith("EXACT_RELEASE")
            and scprint2["scprint2_medium_checkpoint_search"]["artifact_exists"] is False,
            "scPRINT-2 runtime or medium identity differs")
    require(scprint2_exposure["checkpoint_findings"]["scprint2_small_v2"]["exposure_state"]
            == "unknown", "scPRINT-2 exposure differs")
    scprint2_embed = load_json(root / "executions/model-work-302-21075963/embeddings/receipt.json")
    scprint2_eval = load_json(
        root / "executions/model-evaluation-302-21076583/validation/common.stdout.json"
    )
    require(scprint2_embed.get("rows") == 1000
            and scprint2_embed.get("embedding_width") == 424
            and scprint2_embed.get("neighbor_tensors_used") is False
            and scprint2_embed.get("sealed_outcomes_read") is False,
            "scPRINT-2 outcome-blind embedding evidence differs")
    best = scprint2_eval["head_metrics"]["two_layer_mlp"][
        "donor_class_balanced_macro_f1"
    ]
    delta = scprint2_eval["comparisons"]["two_layer_mlp"][
        "delta_donor_class_balanced_macro_f1"
    ]
    require(best == census_models["scprint2_small_v2"]["completed_best_common_macro_f1"]
            and delta == census_models["scprint2_small_v2"]["delta_macro_f1"]
            and delta < 0.0, "scPRINT-2 nonpromotion evidence differs")

    regformer = authorities["config/artifacts/models/regformer/checkpoints.json"]
    regformer_exposure = authorities[
        "config/artifacts/models/regformer/exposure_audit.json"
    ]
    regformer_inspection = evidence[
        "executions/regformer-weights-only-inspection-21066098/weights_only_inspection.json"
    ]
    require(regformer["canonical_checkpoint"]["sha256"] == "UNRESOLVED"
            and regformer_inspection.get("sha256")
            == census_models["regformer"]["observed_checkpoint_sha256"]
            and regformer_inspection.get("tensor_count") == 139,
            "RegFormer observed-versus-registered checkpoint state differs")
    require(regformer_exposure["checkpoint_findings"]["regformer"]["exposure_state"]
            == "unknown", "RegFormer exposure differs")

    scbert = authorities["config/artifacts/models/scbert/checkpoints.json"]
    require(scbert["checkpoint"]["sha256"] == "UNRESOLVED"
            and scbert["checkpoint"]["terms"] == "UNDECLARED"
            and scbert["checkpoint"]["download_status"].startswith("authentication_required"),
            "scBERT checkpoint gate differs")

    baseline = authorities[
        "config/artifacts/models/cell_bulk_classical_baselines/checkpoints.json"
    ]["models"]["celltypist_baseline"]
    require(baseline.get("pretrained_weights") is False
            and baseline.get("task_status") == "registered_mandatory_baseline_not_runnable"
            and baseline.get("runtime_status").startswith("blocked_celltypist_1.7.1_absent"),
            "CellTypist baseline admission differs")

    return {
        "schema_version": "masld-bench-remaining-cell-specialist-audit-v1",
        "status": "pass_fail_closed",
        "census_sha256": CENSUS_SHA256,
        "authority_sha256": AUTHORITY_SHA256,
        "config_sha256": CONFIG_SHA256,
        "evidence_sha256": EVIDENCE_SHA256,
        "artifact_manifest_sha256": ARTIFACT_SHA256,
        "models_reconciled": sorted(expected),
        "next_executable_50k_model": "scimilarity_v1_1",
        "next_action": "outcome_blind_raw_embedding_extraction_only",
        "submission_authorized": False,
        "models_trained": False,
        "models_scored": False,
        "sealed_sources_accessed": [],
        "sealed_champion_created": False,
        "deep_verified_artifacts": sorted(DEEP_VERIFY),
        "legacy_checksum_verified_artifacts": sorted(LEGACY_CHECKSUM_VERIFY),
        "manifest_identity_verified_artifacts": sorted(set(ARTIFACT_SHA256) - DEEP_VERIFY),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve(strict=True)
    if args.output.exists():
        raise CellSpecialistAuditError("output already exists")
    result = audit(root)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "cell_specialist_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
