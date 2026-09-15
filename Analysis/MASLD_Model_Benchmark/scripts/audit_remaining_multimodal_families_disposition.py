#!/usr/bin/env python3
"""Audit the fail-closed remaining multimodal-family disposition."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import tomllib
from typing import Any


SCHEMA_VERSION = "masld-bench-remaining-multimodal-families-disposition-v1"
BLOCKED_STATUS = "FAIL_CLOSED_TOPOLOGY_AND_OUTPUT_CLASS_RECONCILED_NO_NEW_EXECUTION"
MODEL_IDS = (
    "midas",
    "scmomat",
    "stabmap",
    "mira",
    "scmvp",
    "babel",
    "scbutterfly",
    "scjoint",
    "peakvi",
    "epiagent",
    "scooby_onek1k",
    "scooby_epicardioids",
    "scooby_neurips",
)
SOURCE_LINE = re.compile(r"^([0-9a-f]{64})  (/.+)$")


class RemainingMultimodalDispositionError(ValueError):
    """Raised when an authority, receipt, topology, or closed check drifts."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _resolve_member(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise RemainingMultimodalDispositionError("authority member path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise RemainingMultimodalDispositionError(
            f"authority member is a symlink: {relative}"
        )
    member = candidate.resolve(strict=True)
    try:
        member.relative_to(root)
    except ValueError as error:
        raise RemainingMultimodalDispositionError(
            f"authority member escapes benchmark root: {relative}"
        ) from error
    if not member.is_file():
        raise RemainingMultimodalDispositionError(
            f"authority member is not a file: {relative}"
        )
    return member


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RemainingMultimodalDispositionError(
            f"JSON authority is not an object: {path}"
        )
    return value


def _load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    if not isinstance(value, dict):
        raise RemainingMultimodalDispositionError(
            f"TOML authority is not an object: {path}"
        )
    return value


def _verify_source_inventory(
    root: Path, manifest_path: Path, manifest: dict[str, Any]
) -> dict[str, tuple[str, str]]:
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise RemainingMultimodalDispositionError(
            f"artifact inventory differs: {manifest_path}"
        )
    records = [
        record
        for record in artifacts
        if isinstance(record, dict) and record.get("path") == "source.sha256"
    ]
    if len(records) != 1 or not isinstance(records[0].get("sha256"), str):
        raise RemainingMultimodalDispositionError(
            f"source inventory is absent: {manifest_path}"
        )
    source_path = manifest_path.parent / "source.sha256"
    if source_path.is_symlink() or not source_path.is_file():
        raise RemainingMultimodalDispositionError(
            f"source inventory is unavailable: {manifest_path}"
        )
    if _digest(source_path) != records[0]["sha256"]:
        raise RemainingMultimodalDispositionError(
            f"source inventory drifted: {manifest_path}"
        )

    deviations: dict[str, tuple[str, str]] = {}
    for line in source_path.read_text(encoding="utf-8").splitlines():
        match = SOURCE_LINE.fullmatch(line)
        if match is None:
            raise RemainingMultimodalDispositionError(
                f"source inventory line differs: {manifest_path}"
            )
        expected, absolute = match.groups()
        member = Path(absolute)
        if member.is_symlink():
            raise RemainingMultimodalDispositionError(
                f"source member is a symlink: {member}"
            )
        resolved = member.resolve(strict=False)
        try:
            relative = str(resolved.relative_to(root))
        except ValueError as error:
            raise RemainingMultimodalDispositionError(
                f"source member escapes benchmark root: {member}"
            ) from error
        if not resolved.is_file():
            deviations[relative] = (expected, "MISSING")
            continue
        observed = _digest(resolved)
        if observed != expected:
            deviations[relative] = (expected, observed)
    return deviations


def _registry_records(registry: dict[str, Any], wanted: set[str]) -> dict[str, Any]:
    records = {
        record.get("model_id"): record
        for record in registry.get("models", [])
        if isinstance(record, dict) and record.get("model_id") in wanted
    }
    if set(records) != wanted:
        raise RemainingMultimodalDispositionError("registry roster differs")
    return records


def _assert_crosswalk_binding(root: Path, model_id: str) -> None:
    crosswalk_path = _resolve_member(
        root, f"config/artifacts/models/{model_id}/development_crosswalk.json"
    )
    exposure = _load_json(
        _resolve_member(root, f"config/artifacts/models/{model_id}/exposure_audit.json")
    )
    if exposure.get("development_crosswalk_record", {}).get("sha256") != _digest(
        crosswalk_path
    ):
        raise RemainingMultimodalDispositionError(
            f"{model_id} exposure-to-crosswalk binding differs"
        )


def _assert_checkpoint_contracts(root: Path) -> None:
    checkpoints = {
        model_id: _load_json(
            _resolve_member(
                root, f"config/artifacts/models/{model_id}/checkpoints.json"
            )
        )
        for model_id in MODEL_IDS
    }
    midas = checkpoints["midas"]
    if (
        midas.get("code_license") != "MIT"
        or midas.get("implementation", {}).get("pypi_release", {}).get("version")
        != "0.3.0"
        or "RNA-only held rows" not in midas.get("output_contract", {}).get(
            "common_embedding", ""
        )
        or midas.get("safe_loading_contract", {}).get(
            "external_checkpoint_downloaded"
        )
        is not False
    ):
        raise RemainingMultimodalDispositionError("MIDAS authority differs")

    scmomat = checkpoints["scmomat"]
    if (
        scmomat.get("code_license") != "MIT"
        or scmomat.get("implementation", {}).get("implementation_type")
        != "train_locally_transductive"
        or scmomat.get("output_contract", {}).get("profile_status")
        != "latent-only under the frozen RNA-conditioned ATAC registry; reconstructed matrices are not an admitted RNA-only held-donor profile."
    ):
        raise RemainingMultimodalDispositionError("scMoMaT authority differs")

    stabmap = checkpoints["stabmap"]
    if (
        stabmap.get("code_license") != "GPL-2"
        or stabmap.get("implementation", {}).get("implementation_type")
        != "classical_reference_mapping"
        or stabmap.get("output_contract", {}).get("profile_status")
        != "latent-only under the frozen RNA-conditioned ATAC capability registry"
    ):
        raise RemainingMultimodalDispositionError("StabMap authority differs")

    mira = checkpoints["mira"]
    if (
        mira.get("terms_audit", {}).get("result")
        != "blocked_pending_complete_license_text_or_written_confirmation"
        or mira.get("output_contract", {}).get("profile_status") != "latent_only"
        or mira.get("safe_loading_contract", {}).get("external_checkpoint_downloaded")
        is not False
    ):
        raise RemainingMultimodalDispositionError("MIRA authority differs")

    scmvp = checkpoints["scmvp"]
    if (
        scmvp.get("code_license") != "MIT"
        or scmvp.get("output_contract", {}).get("profile_status") != "latent_only"
        or "conditional on observed ATAC" not in scmvp.get("output_contract", {}).get(
            "not_admitted", ""
        )
    ):
        raise RemainingMultimodalDispositionError("scMVP authority differs")

    babel = checkpoints["babel"]
    if (
        babel.get("code_license")
        != "NO_LICENSE_FILE_OR_REPOSITORY_LICENSE_DETECTED_BLOCKED"
        or babel.get("weight_license")
        != "NO_SEPARATE_CHECKPOINT_TERMS_DETECTED_BLOCKED"
        or babel.get("implementation", {})
        .get("checkpoint_identity", {})
        .get("status")
        != "REGISTERED_NOT_DOWNLOADED_UNVERIFIED_SHA256_SIZE_AND_ARCHIVE_CONTENTS"
        or babel.get("safe_loading_contract", {}).get("archive_downloaded")
        is not False
    ):
        raise RemainingMultimodalDispositionError("BABEL authority differs")

    butterfly = checkpoints["scbutterfly"]
    if (
        butterfly.get("code_license") != "MIT"
        or butterfly.get("implementation", {}).get("pypi_release", {}).get("version")
        != "0.0.9"
        or "cell-by-frozen-peak continuous decoder scores"
        not in butterfly.get("output_contract", {}).get("rna_to_atac", "")
        or butterfly.get("safe_loading_contract", {}).get(
            "external_checkpoint_downloaded"
        )
        is not False
    ):
        raise RemainingMultimodalDispositionError("scButterfly authority differs")

    scjoint = checkpoints["scjoint"]
    if (
        scjoint.get("terms_audit", {}).get("result")
        != "blocked_no_detected_license"
        or scjoint.get("output_contract", {}).get("profile_status") != "latent_only"
        or "no peak decoder" not in scjoint.get("output_contract", {}).get(
            "not_admitted", ""
        )
    ):
        raise RemainingMultimodalDispositionError("scJoint authority differs")

    peakvi = checkpoints["peakvi"]
    if (
        peakvi.get("code_license") != "BSD-3-Clause"
        or peakvi.get("implementation", {}).get("pypi_release", {}).get("version")
        != "1.5.0.post1"
        or peakvi.get("safe_loading_contract", {}).get(
            "external_checkpoint_downloaded"
        )
        is not False
    ):
        raise RemainingMultimodalDispositionError("PeakVI authority differs")

    epiagent = checkpoints["epiagent"]
    if (
        epiagent.get("code_license") != "MIT"
        or epiagent.get("weight_license") != "UNRESOLVED"
        or epiagent.get("terms_audit", {}).get("result")
        != "blocked_pending_checkpoint_terms_and_digest"
        or epiagent.get("safe_loading_contract", {}).get("checkpoint_downloaded")
        is not False
    ):
        raise RemainingMultimodalDispositionError("EpiAgent authority differs")

    one = checkpoints["scooby_onek1k"]
    epi = checkpoints["scooby_epicardioids"]
    neu = checkpoints["scooby_neurips"]
    if (
        one.get("safe_loading_contract", {}).get("weight_downloaded") is not False
        or one.get("output_contract", {}).get("profiles") != ["rna_plus", "rna_minus"]
        or "no ATAC output head" not in one.get("output_contract", {}).get(
            "status", ""
        )
        or epi.get("safe_loading_contract", {}).get("weight_downloaded") is not False
        or epi.get("output_contract", {}).get("profiles")
        != ["rna_plus", "rna_minus", "atac_insertion"]
        or "pseudo-matched unpaired" not in epi.get("output_contract", {}).get(
            "warning", ""
        )
        or neu.get("safe_loading_contract", {}).get("weight_downloaded") is not False
        or neu.get("weight_license") != "UNDECLARED"
        or neu.get("output_contract", {}).get("profiles")
        != ["rna_plus", "rna_minus", "atac_insertion"]
    ):
        raise RemainingMultimodalDispositionError("scooby authority differs")


def audit_disposition(root: Path, authority_path: Path) -> dict[str, Any]:
    """Verify metadata-only evidence and keep every new execution check closed."""

    root = root.resolve(strict=True)
    authority_path = authority_path.resolve(strict=True)
    try:
        authority_path.relative_to(root)
    except ValueError as error:
        raise RemainingMultimodalDispositionError(
            "remaining multimodal authority escapes benchmark root"
        ) from error
    authority = _load_json(authority_path)
    if authority.get("schema_version") != SCHEMA_VERSION:
        raise RemainingMultimodalDispositionError("disposition schema differs")
    if authority.get("status") != BLOCKED_STATUS:
        raise RemainingMultimodalDispositionError("disposition is not fail-closed")

    hash_families = {
        "model_authority_hashes": 41,
        "dataset_view_hashes": 6,
        "runtime_manifest_hashes": 5,
        "execution_manifest_hashes": 20,
    }
    manifests: dict[str, dict[str, Any]] = {}
    verified: dict[str, str] = {}
    for family, count in hash_families.items():
        records = authority.get(family)
        if not isinstance(records, dict) or len(records) != count:
            raise RemainingMultimodalDispositionError(f"{family} differs")
        for relative, expected in sorted(records.items()):
            if not isinstance(relative, str) or not isinstance(expected, str):
                raise RemainingMultimodalDispositionError(f"{family} record differs")
            member = _resolve_member(root, relative)
            observed = _digest(member)
            if observed != expected:
                raise RemainingMultimodalDispositionError(
                    f"authority drifted: {relative}"
                )
            verified[relative] = observed
            if family.endswith("manifest_hashes"):
                manifest = _load_json(member)
                if manifest.get("schema_version") != "masld-bench-artifacts-v1":
                    raise RemainingMultimodalDispositionError(
                        f"artifact schema differs: {relative}"
                    )
                manifests[relative] = manifest

    output_classes = authority.get("native_output_classes", {})
    class_keys = (
        "assay_native_rna_to_atac_profile",
        "latent_embedding_or_retrieval_only",
        "observed_atac_latent_only",
        "sequence_to_profile_with_observed_context_not_rna_conditioned_as_released",
    )
    classified = [
        model_id for key in class_keys for model_id in output_classes.get(key, [])
    ]
    if len(classified) != len(set(classified)) or set(classified) != set(MODEL_IDS):
        raise RemainingMultimodalDispositionError("native output partition differs")

    gates = authority.get("task_gates", {})
    if set(gates) != {
        "rna_conditioned_atac",
        "cell_state_mapping",
        "observed_atac_representation",
        "observed_multiome_sequence_context",
    }:
        raise RemainingMultimodalDispositionError("task gate roster differs")
    for gate in gates.values():
        if gate.get("admitted_model_count") != 0:
            raise RemainingMultimodalDispositionError("a model gate opened")
        if gate.get("execution_authorized") is True:
            raise RemainingMultimodalDispositionError("an execution gate opened")
        if gate.get("checkpoint_forward_authorized") is True:
            raise RemainingMultimodalDispositionError("a checkpoint forward opened")
    for gate_id in ("observed_atac_representation", "observed_multiome_sequence_context"):
        if gates[gate_id].get("task_spec_exists") is not False:
            raise RemainingMultimodalDispositionError(
                f"{gate_id} TaskSpec gate differs"
            )

    dispositions = authority.get("model_dispositions", {})
    if set(dispositions) != set(MODEL_IDS) or any(
        record.get("admitted") is not False for record in dispositions.values()
    ):
        raise RemainingMultimodalDispositionError("model disposition opened")

    multimodal_wanted = {
        "midas",
        "scmomat",
        "stabmap",
        "mira",
        "scmvp",
        "babel",
        "scbutterfly",
        "scjoint",
    }
    multimodal = _registry_records(
        _load_toml(_resolve_member(root, "config/models/multimodal_integration.toml")),
        multimodal_wanted,
    )
    expected_multimodal = {
        "midas": ("MIT", "candidate", ["cell_state_mapping", "rna_conditioned_atac"]),
        "scmomat": ("MIT", "deferred", []),
        "stabmap": ("GPL-2.0", "candidate", ["cell_state_mapping"]),
        "mira": (
            "BSD3_declared_complete_license_text_not_detected",
            "blocked_terms",
            ["cell_state_mapping"],
        ),
        "scmvp": ("MIT", "candidate", ["cell_state_mapping"]),
        "babel": (
            "code_NO_LICENSE_DETECTED_weights_NO_TERMS_DETECTED",
            "blocked_terms",
            ["rna_conditioned_atac"],
        ),
        "scbutterfly": (
            "MIT",
            "candidate",
            ["cell_state_mapping", "rna_conditioned_atac"],
        ),
        "scjoint": (
            "blocked_no_detected_license",
            "blocked_terms",
            ["cell_state_mapping"],
        ),
    }
    for model_id, expected in expected_multimodal.items():
        record = multimodal[model_id]
        observed = (
            record.get("license_status"),
            record.get("status"),
            record.get("supported_tasks"),
        )
        if observed != expected or record.get("admission_blocking") is not True:
            raise RemainingMultimodalDispositionError(
                f"{model_id} registry disposition differs"
            )

    chromatin = _registry_records(
        _load_toml(_resolve_member(root, "config/models/chromatin.toml")),
        {"peakvi", "epiagent"},
    )
    expected_chromatin = {
        "peakvi": ("BSD-3-Clause", "deferred", []),
        "epiagent": ("code_MIT_weights_UNDECLARED", "deferred", []),
    }
    for model_id, expected in expected_chromatin.items():
        record = chromatin[model_id]
        observed = (
            record.get("license_status"),
            record.get("status"),
            record.get("supported_tasks"),
        )
        if observed != expected or record.get("admission_blocking") is not True:
            raise RemainingMultimodalDispositionError(
                f"{model_id} registry disposition differs"
            )

    regulatory = _registry_records(
        _load_toml(_resolve_member(root, "config/models/regulatory_sequence.toml")),
        {"scooby_onek1k", "scooby_epicardioids", "scooby_neurips"},
    )
    expected_regulatory = {
        "scooby_onek1k": ("blocked_terms", ["variant_to_regulation"]),
        "scooby_epicardioids": ("blocked_terms", ["variant_to_regulation"]),
        "scooby_neurips": ("restricted_comparator", ["variant_to_regulation"]),
    }
    for model_id, expected in expected_regulatory.items():
        record = regulatory[model_id]
        observed = (record.get("status"), record.get("supported_tasks"))
        if observed != expected or record.get("admission_blocking") is not True:
            raise RemainingMultimodalDispositionError(
                f"{model_id} regulatory disposition differs"
            )

    for model_id in MODEL_IDS:
        _assert_crosswalk_binding(root, model_id)
    _assert_checkpoint_contracts(root)

    view_root = "config/artifacts/dataset-views/gse296875_rna_atac_smoke_1000_v1"
    donor_join = _load_json(_resolve_member(root, f"{view_root}/donor_join.json"))
    topology = _load_json(_resolve_member(root, f"{view_root}/topology.json"))
    reference = _load_json(_resolve_member(root, f"{view_root}/reference.json"))
    qc = _load_json(_resolve_member(root, f"{view_root}/qc.json"))
    exposure = _load_json(_resolve_member(root, f"{view_root}/exposure_audit.json"))
    rights = _load_json(_resolve_member(root, f"{view_root}/rights.json"))
    if (
        donor_join.get("evidence", {}).get("join_key") != "well_id+raw_10x_barcode"
        or donor_join.get("evidence", {}).get("donors") != 39
        or donor_join.get("evidence", {}).get(
            "all_modalities_from_one_donor_stay_in_one_outer_fold"
        )
        is not True
        or topology.get("evidence", {}).get("pairing") != "same_nucleus"
        or topology.get("evidence", {}).get("held_ATAC_available_to_RNA_only_inference")
        is not False
        or reference.get("evidence", {}).get("project_genome") != "GRCh38.p14"
        or reference.get("evidence", {}).get("project_gtf") != "GENCODE_v49"
        or qc.get("evidence", {}).get("selected_nuclei") != 1000
        or qc.get("evidence", {}).get("missing_as_zero") is not False
        or exposure.get("evidence", {}).get("role") != "development_only"
        or exposure.get("evidence", {}).get("sealed_claim_eligible") is not False
        or rights.get("evidence", {}).get("access") != "public_GEO_download"
    ):
        raise RemainingMultimodalDispositionError("dataset topology differs")

    runtime_expectations = {
        "executions/environments/midas-v0.3.0-portable-21066061/ARTIFACTS.json": (
            "midas_v0_3_0_portable_runtime",
            False,
        ),
        "executions/environments/scbutterfly-v0.0.9-portable-21066062/ARTIFACTS.json": (
            "scbutterfly_v0_0_9_portable_runtime",
            False,
        ),
        "executions/environments/scvi-tools-1.5.0.post1-portable-21064579/ARTIFACTS.json": (
            "scvi_tools_1_5_0_post1_portable_runtime",
            True,
        ),
        "executions/environments/stabmap-v1.6.0-portable-21066039/ARTIFACTS.json": (
            "stabmap_v1_6_0_portable_runtime",
            False,
        ),
        "executions/environments/stabmap-v1.6.0-sparse-portable-21066444/ARTIFACTS.json": (
            "stabmap_v1_6_0_sparse_portable_runtime",
            False,
        ),
    }
    for relative, (artifact_class, gpu_pending) in runtime_expectations.items():
        metadata = manifests[relative].get("metadata", {})
        if (
            metadata.get("artifact_class") != artifact_class
            or metadata.get("status") != "passed"
            or metadata.get("outcomes_read") is not False
            or metadata.get("project_data_read") is not False
            or bool(metadata.get("gpu_probe_pending", False)) is not gpu_pending
        ):
            raise RemainingMultimodalDispositionError(
                f"runtime receipt differs: {relative}"
            )

    midas_folds = {
        f"executions/midas-inductive-seed2711-fold{fold}-210{job}/ARTIFACTS.json": fold
        for fold, job in ((0, "68567"), (1, "69966"), (2, "69966"), (3, "69966"), (4, "69966"))
    }
    for relative, fold in midas_folds.items():
        metadata = manifests[relative].get("metadata", {})
        if (
            metadata.get("artifact_class") != "midas_inductive_fit_predict_fold"
            or metadata.get("model_id") != "midas_inductive"
            or metadata.get("fold") != fold
            or metadata.get("seed") != 2711
            or metadata.get("held_atac_input_exposed") is not False
            or metadata.get("query_rna_used_for_training") is not False
            or metadata.get("sealed_inference_eligible") is not False
            or metadata.get("champion_claim_allowed") is not False
            or metadata.get("status") != "passed"
        ):
            raise RemainingMultimodalDispositionError(
                f"MIDAS fold receipt differs: {relative}"
            )
    aggregate_path = (
        "executions/midas-inductive-seed2711-aggregated-21070221/ARTIFACTS.json"
    )
    aggregate = manifests[aggregate_path].get("metadata", {})
    expected_fold_hashes = {
        str(fold): authority["execution_manifest_hashes"][relative]
        for relative, fold in midas_folds.items()
    }
    if (
        aggregate.get("fold_artifacts_sha256") != expected_fold_hashes
        or aggregate.get("prediction_frozen_before_outcomes") is not True
        or aggregate.get("held_atac_input_exposed") is not False
        or aggregate.get("query_rna_used_for_training") is not False
        or aggregate.get("sealed_inference_eligible") is not False
    ):
        raise RemainingMultimodalDispositionError("MIDAS aggregate differs")
    midas_eval = manifests[
        "executions/midas-inductive-seed2711-evaluation-21070222/ARTIFACTS.json"
    ].get("metadata", {})
    if (
        midas_eval.get("prediction_artifacts_sha256")
        != authority["execution_manifest_hashes"][aggregate_path]
        or midas_eval.get("rna_conditioned_atac_eligible") is not True
        or midas_eval.get("held_atac_exposed_to_model") is not False
        or midas_eval.get("prediction_frozen_before_outcomes") is not True
        or midas_eval.get("sealed_inference_eligible") is not False
        or midas_eval.get("seed") != 2711
    ):
        raise RemainingMultimodalDispositionError("MIDAS evaluation differs")

    for relative in (
        "executions/midas-official-trainer-probe-21066691/ARTIFACTS.json",
        "executions/midas-synthetic-training-probe-21066448/ARTIFACTS.json",
    ):
        metadata = manifests[relative].get("metadata", {})
        if (
            metadata.get("model_id") != "midas"
            or metadata.get("synthetic_only") is not True
            or metadata.get("project_data_read") is not False
            or metadata.get("outcomes_read") is not False
            or metadata.get("champion_claim_allowed") is not False
            or metadata.get("status") != "passed"
        ):
            raise RemainingMultimodalDispositionError(
                f"MIDAS probe receipt differs: {relative}"
            )

    butterfly_fit_path = (
        "executions/scbutterfly-b-safe-fit-predict-21066188/ARTIFACTS.json"
    )
    butterfly_fit = manifests[butterfly_fit_path].get("metadata", {})
    butterfly_eval = manifests[
        "executions/scbutterfly-b-safe-evaluation-21066294/ARTIFACTS.json"
    ].get("metadata", {})
    if (
        butterfly_fit.get("model_id") != "scbutterfly_b_safe"
        or butterfly_fit.get("held_atac_input_exposed") is not False
        or butterfly_fit.get("high_level_combined_preprocessing_executed") is not False
        or butterfly_fit.get("inverse_tfidf_executed") is not False
        or butterfly_fit.get("query_atac_state") != "structurally_missing"
        or butterfly_eval.get("prediction_artifacts_sha256")
        != authority["execution_manifest_hashes"][butterfly_fit_path]
        or butterfly_eval.get("rna_conditioned_atac_eligible") is not True
        or butterfly_eval.get("held_atac_exposed_to_model") is not False
        or butterfly_eval.get("sealed_inference_eligible") is not False
    ):
        raise RemainingMultimodalDispositionError("scButterfly evidence differs")

    peak_fit_path = "executions/peakvi-training-context-smoke-21065383/ARTIFACTS.json"
    peak_fit = manifests[peak_fit_path].get("metadata", {})
    peak_eval = manifests[
        "executions/peakvi-training-context-evaluation-21065488/ARTIFACTS.json"
    ].get("metadata", {})
    if (
        peak_fit.get("model_id") != "peakvi_training_context"
        or peak_fit.get("held_atac_available_to_adapter") is not False
        or peak_fit.get("rna_conditioned_atac_eligible") is not False
        or peak_fit.get("sealed_inference_eligible") is not False
        or peak_eval.get("prediction_artifacts_sha256")
        != authority["execution_manifest_hashes"][peak_fit_path]
        or peak_eval.get("rna_conditioned_atac_eligible") is not False
        or peak_eval.get("held_atac_exposed_to_model") is not False
        or peak_eval.get("sealed_inference_eligible") is not False
    ):
        raise RemainingMultimodalDispositionError("PeakVI context evidence differs")

    stab_evaluations = {
        "4817": "executions/stabmap-retrieval-evaluation-21067154/ARTIFACTS.json",
        "5839": "executions/stabmap-retrieval-seed5839-evaluation-21067341/ARTIFACTS.json",
        "6857": "executions/stabmap-retrieval-seed6857-evaluation-21067568/ARTIFACTS.json",
    }
    for relative in stab_evaluations.values():
        metadata = manifests[relative].get("metadata", {})
        if (
            metadata.get("model_id") != "stabmap"
            or metadata.get("hidden_pair_map_available_to_model") is not False
            or metadata.get("cells_used_as_biological_replicates") is not False
            or metadata.get("sealed_rna_conditioned_atac_eligible") is not False
            or metadata.get("champion_claim_allowed") is not False
            or metadata.get("status") != "passed"
        ):
            raise RemainingMultimodalDispositionError(
                f"StabMap evaluation differs: {relative}"
            )
    stability = manifests[
        "executions/stabmap-retrieval-seed-stability-21067940/ARTIFACTS.json"
    ].get("metadata", {})
    expected_evaluations = {
        seed: authority["execution_manifest_hashes"][relative]
        for seed, relative in stab_evaluations.items()
    }
    if (
        stability.get("evaluation_artifacts_sha256") != expected_evaluations
        or stability.get("hidden_pair_map_available_to_models") is not False
        or stability.get("seeds_are_biological_replicates") is not False
        or stability.get("sealed_inference_eligible") is not False
    ):
        raise RemainingMultimodalDispositionError("StabMap stability differs")

    cross = manifests["executions/scooby-cross-review-21064420/ARTIFACTS.json"].get(
        "metadata", {}
    )
    archive = manifests[
        "executions/scooby-released-archive-inventory-21064432/ARTIFACTS.json"
    ].get("metadata", {})
    headers = manifests["executions/scooby-bounded-headers-21064492/ARTIFACTS.json"].get(
        "metadata", {}
    )
    if (
        cross.get("checkpoint_forward_allowed") != []
        or cross.get("rna_conditioned_atac_models") != []
        or cross.get("sealed_inference_allowed") != []
        or archive.get("archive_payload_downloaded") is not False
        or archive.get("checkpoint_bytes_downloaded") is not False
        or headers.get("checkpoint_payload_downloaded") is not False
        or headers.get("checkpoint_deserialized") is not False
        or headers.get("model_forward_executed") is not False
        or headers.get("checkpoint_forward_allowed") != []
        or headers.get("rna_conditioned_atac_models") != []
    ):
        raise RemainingMultimodalDispositionError("scooby bounded evidence differs")

    source_deviations: dict[str, tuple[str, str]] = {}
    for relative, manifest in manifests.items():
        source_deviations.update(
            _verify_source_inventory(root, _resolve_member(root, relative), manifest)
        )
    expected_source_deviations = {
        "executions/environments/stabmap-v1.6.0-sparse-portable-21066444.staging/source/sparseMatrixStats_1.24.0.tar.gz": (
            "4ca8dc253f0630eb4d85a8ba850aa3f3cbb885f3493971a435241d464cfa2564",
            "MISSING",
        ),
        "scripts/scbutterfly_fit_predict_smoke.py": (
            "0fe6a5fd22c16e458489e7853d11b510887b5c7715807f8a349d1d3459229c7c",
            "c61a59a00748d4b7e35d4d40f24a37f0bda654b4a33a7593f4b41241dad523df",
        ),
        "tests/unit/test_scbutterfly_fit_predict_smoke.py": (
            "6bf55dfb86abe775b3838c797ea872d79b5e48f669aed242084e3d7027dd4e7d",
            "eb8715af7d26635bd9d04a147bf2a3c097a6f61a569c03ea8faa1e214f443c18",
        ),
    }
    if source_deviations != expected_source_deviations:
        raise RemainingMultimodalDispositionError(
            f"historical source deviations differ: {source_deviations}"
        )

    drift = authority.get("known_source_and_registry_drift", [])
    if not isinstance(drift, list) or len(drift) != 6:
        raise RemainingMultimodalDispositionError("known drift roster differs")
    drift_by_path = {
        record.get("path"): record for record in drift if isinstance(record, dict)
    }
    if set(drift_by_path) != {
        "executions/environments/stabmap-v1.6.0-sparse-portable-21066444/source.sha256",
        "scripts/scbutterfly_fit_predict_smoke.py",
        "tests/unit/test_scbutterfly_fit_predict_smoke.py",
        "config/evaluation/family_native_tournament.toml",
        "config/evaluation/rna_conditioned_atac_capabilities.toml",
        "config/models/regulatory_sequence.toml",
    }:
        raise RemainingMultimodalDispositionError("known drift paths differ")
    for relative in (
        "scripts/scbutterfly_fit_predict_smoke.py",
        "tests/unit/test_scbutterfly_fit_predict_smoke.py",
        "config/evaluation/family_native_tournament.toml",
        "config/evaluation/rna_conditioned_atac_capabilities.toml",
        "config/models/regulatory_sequence.toml",
    ):
        if _digest(_resolve_member(root, relative)) != drift_by_path[relative].get(
            "observed_sha256"
        ):
            raise RemainingMultimodalDispositionError(
                f"recorded working-tree drift changed: {relative}"
            )
    sparse_member = root / drift_by_path[
        "executions/environments/stabmap-v1.6.0-sparse-portable-21066444/source.sha256"
    ]["expected_member"]
    if sparse_member.exists():
        raise RemainingMultimodalDispositionError(
            "StabMap sparse source member now exists and requires a new audit"
        )

    execution = authority.get("execution_disposition", {})
    false_gates = (
        "dependency_install_authorized",
        "runtime_build_authorized",
        "checkpoint_open_authorized",
        "data_matrix_open_authorized",
        "training_authorized",
        "prediction_authorized",
        "development_scoring_authorized",
        "sealed_scoring_authorized",
        "cpu_campaign_authorized",
        "gpu_campaign_authorized",
        "sbatch_submission_authorized",
    )
    if execution.get("metadata_audit_authorized") is not True or any(
        execution.get(key) is not False for key in false_gates
    ):
        raise RemainingMultimodalDispositionError("execution disposition opened")

    return {
        "status": "pass_metadata_audit_all_new_execution_gates_closed",
        "verified_metadata_files": len(verified),
        "model_identity_count": len(MODEL_IDS),
        "runtime_manifest_count": len(authority["runtime_manifest_hashes"]),
        "execution_manifest_count": len(authority["execution_manifest_hashes"]),
        "historical_source_deviation_count": len(source_deviations),
        "shared_registry_drift_count": 3,
        "admitted_model_count": 0,
        "training_authorized": False,
        "prediction_authorized": False,
        "scoring_authorized": False,
        "sbatch_submission_authorized": False,
        "sealed_outcomes_opened": False,
        "development_outcomes_opened": False,
        "model_payloads_opened": False,
        "data_matrices_opened": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Benchmark repository root.",
    )
    parser.add_argument(
        "--authority",
        type=Path,
        default=None,
        help="Disposition JSON. Defaults below --root.",
    )
    args = parser.parse_args()
    authority = args.authority or (
        args.root
        / "config/artifacts/models/remaining_multimodal_families/disposition_20260824.json"
    )
    receipt = audit_disposition(args.root, authority)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
