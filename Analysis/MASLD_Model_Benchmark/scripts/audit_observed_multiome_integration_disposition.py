#!/usr/bin/env python3
"""Audit the fail-closed observed-multiome integration disposition."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import tomllib
from typing import Any


SCHEMA_VERSION = "masld-bench-observed-multiome-integration-disposition-v2"
BLOCKED_STATUS = "FAIL_CLOSED_TASKSPEC_EXECUTION_BLOCKED_AND_NO_ADMITTED_MODEL_LANE"
CANONICAL_TASK_SPEC = "config/tasks/observed_multiome.toml"
STANDALONE_TASK_SPEC = "config/evaluation/observed_multiome_task.toml"
STANDALONE_DEVELOPMENT_GATE = "config/evaluation/observed_multiome_development_gate.toml"
MODEL_IDS = ("multivi", "scglue", "mofaplus", "seurat_wnn")
SOURCE_LINE = re.compile(r"^([0-9a-f]{64})  (/.+)$")


class ObservedMultiomeDispositionError(ValueError):
    """Raised when an authority, execution receipt, or closed check drifts."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _resolve_member(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ObservedMultiomeDispositionError("authority member path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ObservedMultiomeDispositionError(f"authority member is a symlink: {relative}")
    member = candidate.resolve(strict=True)
    try:
        member.relative_to(root)
    except ValueError as error:
        raise ObservedMultiomeDispositionError(
            f"authority member escapes benchmark root: {relative}"
        ) from error
    if not member.is_file():
        raise ObservedMultiomeDispositionError(
            f"authority member is not a file: {relative}"
        )
    return member


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ObservedMultiomeDispositionError(f"JSON authority is not an object: {path}")
    return value


def _load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    if not isinstance(value, dict):
        raise ObservedMultiomeDispositionError(f"TOML authority is not an object: {path}")
    return value


def _verify_source_inventory(
    root: Path, manifest_path: Path, manifest: dict[str, Any]
) -> dict[str, tuple[str, str]]:
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise ObservedMultiomeDispositionError("artifact inventory differs")
    source_records = [
        record
        for record in artifacts
        if isinstance(record, dict) and record.get("path") == "source.sha256"
    ]
    if len(source_records) != 1 or not isinstance(
        source_records[0].get("sha256"), str
    ):
        raise ObservedMultiomeDispositionError(
            f"source inventory is absent: {manifest_path}"
        )
    record = source_records[0]
    source_path = manifest_path.parent / "source.sha256"
    if source_path.is_symlink() or not source_path.is_file():
        raise ObservedMultiomeDispositionError(
            f"source inventory is unavailable: {manifest_path}"
        )
    if _digest(source_path) != record["sha256"]:
        raise ObservedMultiomeDispositionError(
            f"source inventory drifted: {manifest_path}"
        )
    mismatches: dict[str, tuple[str, str]] = {}
    for line in source_path.read_text(encoding="utf-8").splitlines():
        match = SOURCE_LINE.fullmatch(line)
        if match is None:
            raise ObservedMultiomeDispositionError(
                f"source inventory line differs: {manifest_path}"
            )
        expected, absolute = match.groups()
        member = Path(absolute)
        if member.is_symlink():
            raise ObservedMultiomeDispositionError(
                f"source member is a symlink: {member}"
            )
        resolved = member.resolve(strict=True)
        try:
            relative = str(resolved.relative_to(root))
        except ValueError as error:
            raise ObservedMultiomeDispositionError(
                f"source member escapes benchmark root: {member}"
            ) from error
        observed = _digest(resolved)
        if observed != expected:
            mismatches[relative] = (expected, observed)
    return mismatches


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
        raise ObservedMultiomeDispositionError(
            f"{model_id} exposure-to-crosswalk binding differs"
        )


def audit_disposition(root: Path, authority_path: Path) -> dict[str, Any]:
    """Verify metadata-only evidence and keep every execution check closed."""

    root = root.resolve(strict=True)
    authority_path = authority_path.resolve(strict=True)
    try:
        authority_path.relative_to(root)
    except ValueError as error:
        raise ObservedMultiomeDispositionError(
            "observed-multiome authority escapes benchmark root"
        ) from error
    authority = _load_json(authority_path)
    if authority.get("schema_version") != SCHEMA_VERSION:
        raise ObservedMultiomeDispositionError("disposition schema differs")
    if authority.get("status") != BLOCKED_STATUS:
        raise ObservedMultiomeDispositionError("disposition is not fail-closed")

    hash_families = {
        "model_authority_hashes": 13,
        "dataset_view_hashes": 6,
        "runtime_manifest_hashes": 4,
        "execution_manifest_hashes": 11,
    }
    loaded_manifests: dict[str, dict[str, Any]] = {}
    verified: dict[str, str] = {}
    for family, count in hash_families.items():
        records = authority.get(family)
        if not isinstance(records, dict) or len(records) != count:
            raise ObservedMultiomeDispositionError(f"{family} differs")
        for relative, expected in sorted(records.items()):
            if not isinstance(relative, str) or not isinstance(expected, str):
                raise ObservedMultiomeDispositionError(f"{family} record differs")
            member = _resolve_member(root, relative)
            observed = _digest(member)
            if observed != expected:
                raise ObservedMultiomeDispositionError(f"authority drifted: {relative}")
            verified[relative] = observed
            if family.endswith("manifest_hashes"):
                manifest = _load_json(member)
                if manifest.get("schema_version") != "masld-bench-artifacts-v1":
                    raise ObservedMultiomeDispositionError(
                        f"artifact schema differs: {relative}"
                    )
                loaded_manifests[relative] = manifest

    task_gate = authority.get("task_gate", {})
    if (
        task_gate.get("task_id") != "observed_multiome"
        or task_gate.get("registry_status_observed") != "planned_registration"
        or task_gate.get("canonical_task_spec_exists") is not False
        or task_gate.get("standalone_task_spec_exists") is not True
        or task_gate.get("external_seal_bound") is not False
        or task_gate.get("model_admission_complete") is not False
        or task_gate.get("execution_authorized") is not False
    ):
        raise ObservedMultiomeDispositionError("observed-multiome task gate opened")
    # A standalone TaskSpec was written to config/evaluation on 2026-08-25, after
    # this record.  Its existence does not open the task: it is deliberately kept
    # outside config/tasks so it cannot mutate the frozen global census, and it
    # holds execution closed on its own terms.  So the check is no longer "no
    # TaskSpec exists" but the stricter pair: the canonical registration path is
    # still absent, AND the standalone spec's own checks are still shut.  If
    # either changes, this record needs a new audit.
    expected_paths = {CANONICAL_TASK_SPEC: False, STANDALONE_TASK_SPEC: True}
    recorded_paths = task_gate.get("task_spec_paths")
    if recorded_paths != expected_paths:
        raise ObservedMultiomeDispositionError("TaskSpec path check differs")
    for relative, expected_exists in expected_paths.items():
        if (root / relative).exists() is not expected_exists:
            raise ObservedMultiomeDispositionError(
                f"observed-multiome TaskSpec presence changed at {relative} "
                "and requires a new audit"
            )

    standalone = task_gate.get("standalone_task_spec", {})
    if (
        standalone.get("path") != STANDALONE_TASK_SPEC
        or standalone.get("sha256") != _digest(_resolve_member(root, STANDALONE_TASK_SPEC))
        or standalone.get("registered_in_config_tasks") is not False
        or standalone.get("datasets_sealed") != []
    ):
        raise ObservedMultiomeDispositionError("standalone TaskSpec binding differs")
    spec = _load_toml(_resolve_member(root, STANDALONE_TASK_SPEC))
    if (
        spec.get("task_id") != "observed_multiome"
        or spec.get("status") != standalone.get("document_status")
        or spec.get("datasets_sealed") != []
    ):
        raise ObservedMultiomeDispositionError("standalone TaskSpec drifted")

    recorded_gate = task_gate.get("standalone_development_gate", {})
    if recorded_gate.get("path") != STANDALONE_DEVELOPMENT_GATE or recorded_gate.get(
        "sha256"
    ) != _digest(_resolve_member(root, STANDALONE_DEVELOPMENT_GATE)):
        raise ObservedMultiomeDispositionError("development gate binding differs")
    development_gate = _load_toml(_resolve_member(root, STANDALONE_DEVELOPMENT_GATE))
    if (
        development_gate.get("global_census_revision_required_before_execution") is not True
        or development_gate.get("model_promotion_allowed") is not False
        or development_gate.get("champion_claim_allowed") is not False
        or development_gate.get("selection_lock_allowed") is not False
        or development_gate.get("external_seal_bound") is not False
    ):
        raise ObservedMultiomeDispositionError(
            "the standalone observed-multiome gates opened and require a new audit"
        )

    registry = _load_toml(_resolve_member(root, "config/models/multimodal_integration.toml"))
    model_records = {
        record.get("model_id"): record
        for record in registry.get("models", [])
        if isinstance(record, dict) and record.get("model_id") in MODEL_IDS
    }
    if set(model_records) != set(MODEL_IDS):
        raise ObservedMultiomeDispositionError("multimodal registry roster differs")
    expected_registry = {
        "multivi": (
            "scverse/scvi-tools@56520c713eb1d2b245c72a0f60bc393b74198c91+PyPI:scvi-tools==1.5.0.post1",
            "BSD-3-Clause",
            "candidate",
            ["cell_state_mapping", "rna_conditioned_atac"],
        ),
        "scglue": (
            "gao-lab/GLUE@091bdbbdb3dd6d652ed4677200b034ee3c85dae1+tag:v0.4.1",
            "MIT",
            "candidate",
            ["cell_state_mapping"],
        ),
        "mofaplus": (
            "bioFAM/mofapy2@68c5bee640a643dc4fbfe37c7a0c8001cb318358+PyPI:mofapy2==0.7.4",
            "LGPL-3.0",
            "deferred",
            [],
        ),
        "seurat_wnn": (
            "satijalab/seurat@4c0f2dc16fad4e8f0d7b5c98321d8bcb18caa13a+CRAN:Seurat==5.5.1",
            "MIT",
            "deferred",
            [],
        ),
    }
    for model_id, expected in expected_registry.items():
        record = model_records[model_id]
        observed = (
            record.get("checkpoint_revision"),
            record.get("license_status"),
            record.get("status"),
            record.get("supported_tasks"),
        )
        if observed != expected or record.get("admission_blocking") is not True:
            raise ObservedMultiomeDispositionError(
                f"{model_id} registry disposition differs"
            )

    for model_id in MODEL_IDS:
        _assert_crosswalk_binding(root, model_id)
    multivi = _load_json(
        _resolve_member(root, "config/artifacts/models/multivi/checkpoints.json")
    )
    scglue = _load_json(
        _resolve_member(root, "config/artifacts/models/scglue/checkpoints.json")
    )
    mofaplus = _load_json(
        _resolve_member(root, "config/artifacts/models/mofaplus/checkpoints.json")
    )
    seurat = _load_json(
        _resolve_member(root, "config/artifacts/models/seurat_wnn/checkpoints.json")
    )
    if (
        multivi.get("code_license") != "BSD-3-Clause"
        or multivi.get("implementation", {}).get("pypi_release", {}).get("version")
        != "1.5.0.post1"
        or multivi.get("missingness_contract", {}).get("admission_status")
        != "BLOCKED_PENDING_EXPLICIT_MASK_ADAPTER_OR_APPROVED_CONTRACT_EXCEPTION"
        or multivi.get("safe_loading_contract", {}).get("external_checkpoint_downloaded")
        is not False
        or multivi.get("weight_license") != "NOT_APPLICABLE_PROJECT_TRAINED"
    ):
        raise ObservedMultiomeDispositionError("MultiVI authority differs")
    if (
        scglue.get("code_license") != "MIT"
        or scglue.get("implementation", {}).get("repository_tag") != "v0.4.1"
        or scglue.get("output_contract", {}).get("current_profile_status")
        != "latent-only under the current RNA-conditioned ATAC capability registry"
        or "dill" not in scglue.get("safe_loading_contract", {}).get("upstream_risk", "")
        or scglue.get("weight_license") != "NOT_APPLICABLE_PROJECT_TRAINED"
    ):
        raise ObservedMultiomeDispositionError("scGLUE authority differs")
    if (
        mofaplus.get("code_license") != "LGPL-3.0"
        or mofaplus.get("implementation", {}).get("pypi_release", {}).get("version")
        != "0.7.4"
        or mofaplus.get("output_contract", {}).get("profile_status") != "latent_only"
        or "no amortized encoder" not in mofaplus.get("architecture_contract", {}).get(
            "inductive_limitation", ""
        )
        or mofaplus.get("weight_license")
        != "NOT_APPLICABLE_PROJECT_FITTED_HDF5"
    ):
        raise ObservedMultiomeDispositionError("MOFA+ authority differs")
    if (
        seurat.get("admission_disposition", {}).get("status")
        != "deferred_latent_only_no_registered_task"
        or seurat.get("registered_identity", {}).get("registry_supported_tasks") != []
        or seurat.get("checkpoint_contract", {}).get("trained_artifact_identity")
        != "NOT_YET_CREATED"
        or seurat.get("runtime_contract", {}).get("status")
        != "UNRESOLVED_exact_R_and_compiled_dependency_lock_not_created"
        or seurat.get("terms_audit", {}).get("code")
        != "The pinned repository and CRAN metadata declare MIT plus the package LICENSE file."
    ):
        raise ObservedMultiomeDispositionError("Seurat WNN authority differs")

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
        or reference.get("evidence", {}).get("project_gtf_sha256")
        != "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
        or qc.get("evidence", {}).get("selected_nuclei") != 1000
        or qc.get("evidence", {}).get("missing_as_zero") is not False
        or exposure.get("evidence", {}).get("role") != "development_only"
        or exposure.get("evidence", {}).get("sealed_claim_eligible") is not False
        or rights.get("evidence", {}).get("redistribution")
        != "source_reference_only"
    ):
        raise ObservedMultiomeDispositionError("development fixture authority differs")

    runtime_expected = {
        "executions/environments/scvi-tools-1.5.0.post1-portable-21064579/ARTIFACTS.json": {
            "artifact_class": "scvi_tools_1_5_0_post1_portable_runtime",
            "status": "passed",
            "scvi_tools_version": "1.5.0.post1",
            "gpu_probe_pending": True,
            "outcomes_read": False,
            "project_data_read": False,
        },
        "executions/environments/pytorch-2.6.0-cu124-portable-21064747/ARTIFACTS.json": {
            "artifact_class": "pytorch_2_6_0_cu124_portable_overlay",
            "status": "passed",
            "torch_version": "2.6.0+cu124",
            "gpu_probe_pending": True,
            "outcomes_read": False,
            "project_data_read": False,
        },
        "executions/environments/scglue-v0.4.1-portable-21064673/ARTIFACTS.json": {
            "artifact_class": "scglue_v0_4_1_portable_runtime",
            "status": "passed",
            "scglue_version": "0.4.1",
            "gpu_probe_pending": True,
            "unsafe_dill_load_called": False,
            "outcomes_read": False,
            "project_data_read": False,
        },
        "executions/environments/seurat-5.5.1-portable-21064881/ARTIFACTS.json": {
            "artifact_class": "seurat_5_5_1_portable_runtime",
            "status": "passed",
            "seurat_version": "5.5.1",
            "source_sha256": "9614ef02d3e1010c40be5916a309103a76c4221a667cbc4b312e5126459a5821",
            "model_fit_executed": False,
            "outcomes_read": False,
            "project_data_read": False,
        },
    }
    source_mismatches: dict[str, tuple[str, str]] = {}
    for relative, expected in runtime_expected.items():
        manifest = loaded_manifests[relative]
        metadata = manifest.get("metadata")
        if not isinstance(metadata, dict) or any(
            metadata.get(key) != value for key, value in expected.items()
        ):
            raise ObservedMultiomeDispositionError(f"runtime metadata differs: {relative}")
        source_mismatches.update(
            _verify_source_inventory(root, _resolve_member(root, relative), manifest)
        )

    execution_expected = {
        "executions/multivi-peakvi-no-outcome-probe-21064846/ARTIFACTS.json": {
            "artifact_class": "multivi_peakvi_no_outcome_runtime_probe",
            "status": "passed",
            "champion_claim_allowed": False,
            "multivi_explicit_state_bridge_passed": True,
            "multivi_project_fit_pending": True,
            "outcomes_read": False,
            "project_data_read": False,
        },
        "executions/multivi-smoke-prepared-21065093/ARTIFACTS.json": {
            "artifact_class": "multivi_smoke_prepared_folds",
            "status": "passed",
            "champion_claim_allowed": False,
            "pairing_topology": "same_nucleus",
            "held_atac_exported": False,
            "outcomes_read": False,
        },
        "executions/multivi-smoke-fit-predict-21065203/ARTIFACTS.json": {
            "artifact_class": "multivi_smoke_fit_predict",
            "status": "passed",
            "champion_claim_allowed": False,
            "outer_unit": "donor",
            "query_atac_read": False,
            "query_atac_state": "structurally_missing",
            "outcomes_read": False,
        },
        "executions/multivi-smoke-evaluation-21065278/ARTIFACTS.json": {
            "artifact_class": "multivi_smoke_development_evaluation",
            "status": "passed",
            "champion_claim_allowed": False,
            "held_atac_exposed_to_model": False,
            "prediction_frozen_before_outcomes": True,
            "test_outcomes_read": False,
        },
        "executions/scglue-no-outcome-probe-21064835/ARTIFACTS.json": {
            "artifact_class": "scglue_no_outcome_runtime_probe",
            "status": "passed",
            "champion_claim_allowed": False,
            "latent_integration_only": True,
            "outcomes_read": False,
            "project_data_read": False,
        },
        "executions/scglue-retrieval-prepared-21065511/ARTIFACTS.json": {
            "artifact_class": "scglue_retrieval_prepared_folds",
            "status": "passed",
            "champion_claim_allowed": False,
            "pairing_topology": "same_nucleus",
            "hidden_pair_map_available_to_model": False,
            "retrieval_metrics_calculated": False,
        },
        "executions/scglue-retrieval-fit-predict-21065540/ARTIFACTS.json": {
            "artifact_class": "scglue_retrieval_fit_predict",
            "status": "passed",
            "champion_claim_allowed": False,
            "hidden_pair_map_read": False,
            "retrieval_metrics_calculated": False,
            "rna_to_atac_prediction_executed": False,
            "sealed_rna_conditioned_atac_eligible": False,
        },
        "executions/scglue-retrieval-evaluation-21065607/ARTIFACTS.json": {
            "artifact_class": "scglue_same_nucleus_retrieval_development_evaluation",
            "status": "passed",
            "champion_claim_allowed": False,
            "hidden_pair_map_available_to_models": False,
            "cells_used_as_biological_replicates": False,
            "rna_to_atac_prediction_claim": False,
            "test_outcomes_read": False,
        },
        "executions/seurat-wnn-bridge-prepared-21065772/ARTIFACTS.json": {
            "artifact_class": "seurat_wnn_bridge_prepared_inputs",
            "status": "passed",
            "champion_claim_allowed": False,
            "pairing_topology": "same_nucleus",
            "hidden_pair_map_read": False,
            "rna_to_atac_prediction_executed": False,
            "outcomes_read": False,
        },
        "executions/seurat-wnn-bridge-fit-predict-21065917/ARTIFACTS.json": {
            "artifact_class": "seurat_wnn_bridge_fit_predict",
            "status": "passed",
            "champion_claim_allowed": False,
            "model_id": "seurat_wnn_bridge",
            "hidden_pair_map_read": False,
            "retrieval_metrics_calculated": False,
            "rna_to_atac_prediction_executed": False,
            "outcomes_read": False,
        },
        "executions/seurat-wnn-bridge-evaluation-21065970/ARTIFACTS.json": {
            "artifact_class": "seurat_wnn_bridge_retrieval_development_evaluation",
            "status": "passed",
            "champion_claim_allowed": False,
            "hidden_pair_map_available_to_model": False,
            "cells_used_as_biological_replicates": False,
            "rna_to_atac_prediction_claim": False,
            "test_outcomes_read": False,
        },
    }
    for relative, expected in execution_expected.items():
        manifest = loaded_manifests[relative]
        metadata = manifest.get("metadata")
        if not isinstance(metadata, dict) or any(
            metadata.get(key) != value for key, value in expected.items()
        ):
            raise ObservedMultiomeDispositionError(
                f"execution metadata differs: {relative}"
            )
        source_mismatches.update(
            _verify_source_inventory(root, _resolve_member(root, relative), manifest)
        )

    expected_source_drift = {
        "scripts/evaluate_multivi_smoke.py": (
            "57c0195f55254ad725eea04aaee074c2ac164d5ca038027b7f7fd7cb40939b90",
            "085f973e25a3e6dd83474afe7c8b3bd4d136ae70bee92ea687313c6883cac09f",
        )
    }
    if source_mismatches != expected_source_drift:
        raise ObservedMultiomeDispositionError("execution source-chain drift differs")

    drift = authority.get("known_source_and_registry_drift")
    if not isinstance(drift, list) or len(drift) != 3:
        raise ObservedMultiomeDispositionError("known drift roster differs")
    for finding in drift:
        if not isinstance(finding, dict):
            raise ObservedMultiomeDispositionError("known drift record differs")
        path = finding.get("path")
        observed = finding.get("observed_sha256")
        if not isinstance(path, str) or not isinstance(observed, str):
            raise ObservedMultiomeDispositionError("known drift binding differs")
        if _digest(_resolve_member(root, path)) != observed:
            raise ObservedMultiomeDispositionError(f"known drift changed: {path}")
        embedded = finding.get("embedded_execution_sha256")
        if embedded is None:
            embedded = finding.get("embedded_repository_sha256")
        if embedded == observed:
            raise ObservedMultiomeDispositionError(
                f"known drift unexpectedly closed: {path}"
            )

    if any(root.joinpath("executions").glob("*mofaplus*")):
        raise ObservedMultiomeDispositionError(
            "MOFA+ execution evidence now exists and requires a new audit"
        )
    dispositions = authority.get("model_dispositions", {})
    if set(dispositions) != set(MODEL_IDS) or any(
        dispositions[model_id].get("admitted") is not False for model_id in MODEL_IDS
    ):
        raise ObservedMultiomeDispositionError("model admission gate opened")
    execution = authority.get("execution_disposition", {})
    if execution.get("metadata_audit_authorized") is not True or any(
        execution.get(field) is not False
        for field in (
            "dependency_install_authorized",
            "runtime_build_authorized",
            "checkpoint_or_rds_open_authorized",
            "data_matrix_open_authorized",
            "training_authorized",
            "prediction_authorized",
            "development_scoring_authorized",
            "sealed_scoring_authorized",
            "gpu_sbatch_authorized",
            "observed_multiome_cpu_or_gpu_campaign_authorized",
        )
    ):
        raise ObservedMultiomeDispositionError("execution gate opened")

    return {
        "schema_version": "masld-bench-observed-multiome-disposition-audit-receipt-v1",
        "status": "pass_metadata_audit_all_model_lanes_blocked",
        "authority_sha256": _digest(authority_path),
        "verified_metadata_files": len(verified),
        "runtime_manifest_count": 4,
        "execution_manifest_count": 11,
        "historical_source_drift_count": 1,
        "shared_registry_drift_count": 2,
        "task_spec_bound": False,
        "canonical_task_spec_exists": False,
        "standalone_task_spec_exists": True,
        "standalone_task_spec_execution_authorized": False,
        "external_seal_bound": False,
        "admitted_model_count": 0,
        "training_authorized": False,
        "scoring_authorized": False,
        "gpu_sbatch_authorized": False,
        "sealed_outcomes_opened": False,
        "development_outcomes_opened": False,
        "model_artifacts_opened": False,
        "data_matrices_opened": False,
    }


def main() -> None:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument(
        "--authority",
        type=Path,
        default=default_root
        / "config/artifacts/models/observed_multiome_integration/disposition_20260824.json",
    )
    arguments = parser.parse_args()
    receipt = audit_disposition(arguments.root, arguments.authority)
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
