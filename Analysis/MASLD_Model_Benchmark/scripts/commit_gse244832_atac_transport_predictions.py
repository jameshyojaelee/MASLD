#!/usr/bin/env python3
"""Commit observed-ATAC GSE244832 transport predictions before scoring."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from scripts.commit_gse281367_atac_transport_predictions import (
    COMMIT_SCHEMA,
    OUTPUT_FAMILIES,
    ROTATIONS,
    SCHEMA,
    ATACPredictionCommitError,
    _digest,
    _load_json,
)


DATASET_ID = "gse244832"
AXIS_CLASS = "gse244832_reference_guarded_atac_exchange_axis"
QUERY_CLASS = "gse244832_label_free_query_atac"
REGISTRATION_CLASS = "gse244832_observed_atac_execution_registration"
LAYOUT = "donor_lineage"


def validate_evidence_contract(manifest: Mapping[str, Any]) -> None:
    evidence = manifest.get("evidence_contract")
    if not isinstance(evidence, dict):
        raise ATACPredictionCommitError("prediction evidence contract is absent")
    for field in (
        "condition_labels_consumed",
        "phenotype_values_consumed",
        "rna_assay_consumed",
        "cross_assay_join_consumed",
        "sequence_features_consumed",
        "source_native_reference_bundle_assumed_resolved",
        "development_outcomes_consumed",
        "sealed_data_consumed",
        "target_role_atac_consumed",
        "target_role_contigs_entered_query_transform",
        "missing_as_zero",
    ):
        if evidence.get(field) is not False:
            raise ATACPredictionCommitError(f"prediction firewall opened: {field}")
    context_role, target_role = ROTATIONS.get(str(manifest.get("rotation_id")), (None, None))
    if context_role is None or evidence.get("scored_target_role") != target_role:
        raise ATACPredictionCommitError("rotation or scored target role differs")
    if (
        manifest.get("input_regime") != "observed_atac"
        or manifest.get("prediction_layout") != LAYOUT
        or evidence.get("observed_atac_context_role") != context_role
        or evidence.get("query_atac_consumed") is not True
    ):
        raise ATACPredictionCommitError("GSE244832 permits observed-ATAC donor-lineage views only")
    receptive_field = evidence.get("receptive_field_bp")
    if not (
        receptive_field == "not_applicable_global_observed_atac_encoder"
        or isinstance(receptive_field, int)
        and 1 <= receptive_field <= 524288
    ):
        raise ATACPredictionCommitError("receptive-field declaration differs")
    query_scope = evidence.get("query_transform_scope")
    if query_scope not in {
        "source_training_only",
        "per_donor_context_only",
    }:
        raise ATACPredictionCommitError("query-transform scope differs")
    operation = evidence.get("per_donor_context_operation")
    if (
        query_scope == "source_training_only"
        and operation != "not_applicable"
        or query_scope == "per_donor_context_only"
        and operation != "deterministic_normalization_without_fitted_target_state"
    ):
        raise ATACPredictionCommitError("per-donor query operation differs")


def validate_array_contract(
    values: Any,
    missing: Any,
    unit_states: Any,
    *,
    output_family: str,
    minimum_coverage: float = 0.95,
) -> dict[str, Any]:
    import numpy as np

    score = np.asarray(values)
    states = np.asarray(missing)
    unit_missing = np.asarray(unit_states)
    profile = output_family.endswith("profile")
    expected = (18, 4, 16000, 20) if profile else (18, 4, 16000)
    expected_missing = expected[:-1] if profile else expected
    if (
        tuple(score.shape) != expected
        or tuple(states.shape) != expected_missing
        or unit_missing.shape != (18, 4)
    ):
        raise ATACPredictionCommitError("prediction or missing-state shape differs")
    if (
        score.dtype != np.float32
        or states.dtype != np.uint8
        or unit_missing.dtype != np.uint8
        or np.any(states > 7)
        or not set(np.unique(unit_missing)).issubset({0, 1, 3})
        or int((unit_missing == 0).sum()) != 48
    ):
        raise ATACPredictionCommitError("prediction or missing-state dtype differs")
    for donor in range(18):
        for lineage in range(4):
            fixed_state = int(unit_missing[donor, lineage])
            row_states = states[donor, lineage]
            if fixed_state == 0:
                if np.any(np.isin(row_states, (1, 2, 3))):
                    raise ATACPredictionCommitError("eligible unit uses biological-unit missingness")
            elif not np.all(row_states == fixed_state):
                raise ATACPredictionCommitError("ineligible unit was not masked before prediction")
    observed = states == 0
    expanded = np.broadcast_to(observed[..., None], score.shape) if profile else observed
    if np.any(~np.isfinite(score[expanded])) or np.any(score[expanded] < 0):
        raise ATACPredictionCommitError("observed predictions are not finite and nonnegative")
    if np.any(~np.isnan(score[~expanded])):
        raise ATACPredictionCommitError("missing predictions must be NaN, never zero")
    coverage = observed[unit_missing == 0].reshape(48, 16000).mean(axis=1)
    if np.any(coverage < minimum_coverage):
        raise ATACPredictionCommitError("prediction coverage is below the frozen threshold")
    return {
        "minimum_coverage": float(coverage.min()),
        "maximum_coverage": float(coverage.max()),
        "observed_predictions": int(observed.sum()),
        "missing_predictions": int((~observed).sum()),
        "eligible_units": int((unit_missing == 0).sum()),
        "structurally_missing_units": int((unit_missing == 1).sum()),
        "below_qc_units": int((unit_missing == 3).sum()),
    }


def load_registration(
    root: Path,
    artifacts_sha256: str,
    axis_artifacts_sha256: str,
) -> tuple[dict[str, Any], Any]:
    import numpy as np

    root = root.resolve(strict=True)
    manifest = verify_frozen_tree(root)
    metadata = manifest.get("metadata", {})
    contract = _load_json(root / "execution_contract.json")
    if (
        root.is_symlink()
        or _digest(root / "ARTIFACTS.json") != artifacts_sha256
        or metadata.get("artifact_class") != REGISTRATION_CLASS
        or metadata.get("dataset_id") != DATASET_ID
        or metadata.get("axis_artifacts_sha256") != axis_artifacts_sha256
        or metadata.get("query_lineage_mask_required") is not True
        or metadata.get("biological_outcome_arrays_read") is not False
        or metadata.get("sequence_execution_authorized") is not False
        or contract.get("axis_artifacts_sha256") != axis_artifacts_sha256
        or contract.get("direct_axis_family_table_execution_authorized") is not False
        or contract.get("direct_query_artifact_execution_authorized") is not False
        or contract.get("prediction_commit_required_before_outcome_open") is not True
        or contract.get("sequence_execution_authorized") is not False
        or contract.get("champion_claim_allowed") is not False
    ):
        raise ATACPredictionCommitError("GSE244832 execution registration differs")
    mask_path = (root / str(contract.get("query_lineage_missing_state_path", ""))).resolve(
        strict=True
    )
    mask_path.relative_to(root)
    if (
        mask_path.is_symlink()
        or _digest(mask_path) != contract.get("query_lineage_missing_state_sha256")
    ):
        raise ATACPredictionCommitError("query lineage missingness binding differs")
    unit_states = np.load(mask_path, mmap_mode="r", allow_pickle=False)
    if unit_states.shape != (18, 4) or unit_states.dtype != np.uint8:
        raise ATACPredictionCommitError("query lineage missingness geometry differs")
    return contract, unit_states


def validate_query_authority(
    root: Path,
    artifacts_sha256: str,
    context_role: str,
    target_role: str,
) -> None:
    manifest = verify_frozen_tree(root)
    metadata = manifest.get("metadata", {})
    contract = _load_json(root / "query_contract.json")
    if (
        _digest(root / "ARTIFACTS.json") != artifacts_sha256
        or metadata.get("artifact_class") != QUERY_CLASS
        or metadata.get("dataset_id") != DATASET_ID
        or metadata.get("context_role") != context_role
        or metadata.get("scored_target_role") != target_role
        or metadata.get("condition_values_read") is not False
        or metadata.get("rna_assay_opened") is not False
        or metadata.get("evaluator_outcome_artifact_read") is not False
        or metadata.get("sequence_execution_authorized") is not False
        or contract.get("scored_target_contigs_exported_in_query_values") is not False
        or contract.get("cross_assay_join_inferred") is not False
        or contract.get("metrics_calculated") is not False
    ):
        raise ATACPredictionCommitError("query-ATAC authority differs")


def resolve_model_manifest(
    prediction_root: Path,
    record: Mapping[str, Any],
) -> Path:
    relative = Path(str(record.get("path", "")))
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise ATACPredictionCommitError("model manifest escaped prediction authority")
    path = (prediction_root / relative).resolve(strict=True)
    try:
        path.relative_to(prediction_root)
    except ValueError as error:
        raise ATACPredictionCommitError("model manifest escaped prediction authority") from error
    if (
        path.is_symlink()
        or not path.is_file()
        or _digest(path) != record.get("sha256")
    ):
        raise ATACPredictionCommitError("model manifest binding differs")
    return path


def validate_bundle(
    prediction_root: Path,
    prediction_artifacts_sha256: str,
    axis_root: Path,
    axis_artifacts_sha256: str,
    registration_root: Path,
    registration_artifacts_sha256: str,
) -> dict[str, Any]:
    import numpy as np

    prediction_root = prediction_root.resolve(strict=True)
    axis_root = axis_root.resolve(strict=True)
    prediction_artifacts = verify_frozen_tree(prediction_root)
    axis_artifacts = verify_frozen_tree(axis_root)
    axis_summary = _load_json(axis_root / "exchange_summary.json")
    if _digest(prediction_root / "ARTIFACTS.json") != prediction_artifacts_sha256:
        raise ATACPredictionCommitError("prediction tree hash differs")
    if (
        _digest(axis_root / "ARTIFACTS.json") != axis_artifacts_sha256
        or axis_artifacts.get("metadata", {}).get("artifact_class") != AXIS_CLASS
        or axis_artifacts.get("metadata", {}).get("sequence_execution_authorized") is not False
        or axis_summary.get("non_observed_atac_executable_rows") != 0
    ):
        raise ATACPredictionCommitError("reference-guarded exchange-axis authority differs")
    registration, unit_states = load_registration(
        registration_root,
        registration_artifacts_sha256,
        axis_artifacts_sha256,
    )
    manifest = _load_json(prediction_root / "prediction_bundle.json")
    model_id = str(manifest.get("model_id", ""))
    if (
        manifest.get("schema_version") != SCHEMA
        or manifest.get("input_regime") != "observed_atac"
        or manifest.get("output_family") not in OUTPUT_FAMILIES
        or manifest.get("prediction_layout") != LAYOUT
        or manifest.get("axis_artifacts_sha256") != axis_artifacts_sha256
        or manifest.get("target_dataset_id") != DATASET_ID
        or manifest.get("source_dataset_id") != "gse296875"
        or model_id not in registration.get("allowed_model_ids", ())
        or manifest.get("output_family")
        not in registration.get("allowed_output_families_by_model", {}).get(model_id, ())
        or manifest.get("execution_registration_artifacts_sha256")
        != registration_artifacts_sha256
    ):
        raise ATACPredictionCommitError("prediction manifest identity differs")
    validate_evidence_contract(manifest)
    context_role, target_role = ROTATIONS[str(manifest["rotation_id"])]
    query = manifest.get("query_atac_authority")
    if not isinstance(query, dict):
        raise ATACPredictionCommitError("query-ATAC authority binding is absent")
    query_root = Path(str(query.get("path", ""))).resolve(strict=True)
    registered_query = registration.get("query_authorities", {}).get(context_role, {})
    if (
        query_root.is_symlink()
        or query.get("context_role") != context_role
        or str(query_root) != registered_query.get("path")
        or query.get("artifacts_sha256") != registered_query.get("artifacts_sha256")
    ):
        raise ATACPredictionCommitError("query-ATAC path or role differs")
    validate_query_authority(
        query_root,
        str(query.get("artifacts_sha256", "")),
        context_role,
        target_role,
    )
    model_manifest = manifest.get("model_manifest")
    if not isinstance(model_manifest, dict):
        raise ATACPredictionCommitError("model manifest binding is absent")
    model_path = resolve_model_manifest(prediction_root, model_manifest)
    model_contract = _load_json(model_path)
    evidence = manifest["evidence_contract"]
    if (
        model_contract.get("model_id") != model_id
        or model_contract.get("execution_registration_artifacts_sha256")
        != registration_artifacts_sha256
        or model_contract.get("target_dataset_id") != DATASET_ID
        or model_contract.get("fit_scope") != "gse296875_source_training_only"
        or model_contract.get("sequence_features_consumed") is not False
        or model_contract.get("rna_assay_consumed") is not False
        or model_contract.get("target_dataset_fit_or_adaptation") is not False
        or model_contract.get("per_donor_context_operation")
        != evidence.get("per_donor_context_operation")
        or model_id in {"lsi", "observed_atac_glm", "observed_atac_only"}
        and evidence.get("query_transform_scope") != "source_training_only"
    ):
        raise ATACPredictionCommitError("model-specific execution manifest differs")
    arrays = manifest.get("arrays")
    if not isinstance(arrays, dict) or set(arrays) != {"predictions", "missing_state"}:
        raise ATACPredictionCommitError("prediction array roster differs")
    paths: dict[str, Path] = {}
    for name, expected_dtype in (("predictions", "float32"), ("missing_state", "uint8")):
        record = arrays[name]
        path = (prediction_root / str(record.get("path", ""))).resolve(strict=True)
        path.relative_to(prediction_root)
        if path.is_symlink() or _digest(path) != record.get("sha256") or record.get("dtype") != expected_dtype:
            raise ATACPredictionCommitError(f"{name} array binding differs")
        paths[name] = path
    predictions = np.load(paths["predictions"], mmap_mode="r", allow_pickle=False)
    missing = np.load(paths["missing_state"], mmap_mode="r", allow_pickle=False)
    array_receipt = validate_array_contract(
        predictions,
        missing,
        unit_states,
        output_family=str(manifest["output_family"]),
    )
    return {
        "schema_version": COMMIT_SCHEMA,
        "status": "committed_before_outcome_open",
        "bundle_id": manifest.get("bundle_id"),
        "model_id": manifest.get("model_id"),
        "input_regime": "observed_atac",
        "output_family": manifest.get("output_family"),
        "prediction_layout": LAYOUT,
        "rotation_id": manifest.get("rotation_id"),
        "prediction_artifacts_sha256": prediction_artifacts_sha256,
        "axis_artifacts_sha256": axis_artifacts_sha256,
        "execution_registration_artifacts_sha256": registration_artifacts_sha256,
        "query_lineage_missing_state_sha256": registration[
            "query_lineage_missing_state_sha256"
        ],
        "query_atac_artifacts_sha256": query["artifacts_sha256"],
        "model_manifest_sha256": model_manifest["sha256"],
        "array_validation": array_receipt,
        "prediction_values_read_for_structure_validation": True,
        "condition_labels_read": False,
        "phenotype_values_read": False,
        "rna_assay_read": False,
        "development_outcomes_read": False,
        "metrics_calculated": False,
        "sealed_data_read": False,
        "sequence_execution_authorized": False,
        "thresholds_or_calibration_changed": False,
    }


def commit(
    prediction_root: Path,
    prediction_artifacts_sha256: str,
    axis_root: Path,
    axis_artifacts_sha256: str,
    registration_root: Path,
    registration_artifacts_sha256: str,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise ATACPredictionCommitError("refusing to overwrite prediction commit")
    receipt = validate_bundle(
        prediction_root,
        prediction_artifacts_sha256,
        axis_root,
        axis_artifacts_sha256,
        registration_root,
        registration_artifacts_sha256,
    )
    output.mkdir(parents=True, mode=0o750)
    write_json_exclusive(output / "commit.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse244832_atac_transport_prediction_commit",
            "bundle_id": receipt["bundle_id"],
            "model_id": receipt["model_id"],
            "prediction_artifacts_sha256": prediction_artifacts_sha256,
            "axis_artifacts_sha256": axis_artifacts_sha256,
            "execution_registration_artifacts_sha256": registration_artifacts_sha256,
            "development_outcomes_read": False,
            "metrics_calculated": False,
            "sequence_execution_authorized": False,
            "status": "passed",
        },
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction-root", required=True, type=Path)
    parser.add_argument("--prediction-artifacts-sha256", required=True)
    parser.add_argument("--axis-root", required=True, type=Path)
    parser.add_argument("--axis-artifacts-sha256", required=True)
    parser.add_argument("--registration-root", required=True, type=Path)
    parser.add_argument("--registration-artifacts-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    print(json.dumps(commit(**vars(arguments)), sort_keys=True))


if __name__ == "__main__":
    main()
