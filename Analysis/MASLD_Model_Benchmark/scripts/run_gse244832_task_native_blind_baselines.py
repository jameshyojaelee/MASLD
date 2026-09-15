#!/usr/bin/env python3
"""Fit three source-only baselines, predict GSE244832 blindly, and commit."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import (
    freeze_tree,
    verify_frozen_tree,
    write_json_exclusive,
    write_text_exclusive,
)
from masld_bench.gse244832_task_native_baselines import (
    LINEAGES,
    MODEL_IDS,
    export_state,
    fit,
    predict,
)
from scripts.commit_gse244832_atac_transport_predictions import commit


SCHEMA = "masld-bench-gse244832-task-native-blind-baselines-v1"
ROTATIONS = {
    "valid_context_test_target": ("valid", "test"),
    "test_context_valid_target": ("test", "valid"),
}
ROLE_INDEX = {"valid": 0, "test": 1}


class BlindBaselineRunError(ValueError):
    """Raised when source fitting, blind prediction, or commit differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise BlindBaselineRunError(f"JSON object required: {path}")
    return value


def resolve_tree(
    root: Path,
    record: Mapping[str, Any],
    artifact_class: str,
) -> Path:
    path = (root / str(record.get("path", ""))).resolve(strict=True)
    path.relative_to(root)
    manifest = verify_frozen_tree(path)
    if (
        path.is_symlink()
        or digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256")
        or manifest.get("metadata", {}).get("artifact_class") != artifact_class
    ):
        raise BlindBaselineRunError(f"{artifact_class} authority differs")
    return path


def validate_config(root: Path, config_path: Path) -> tuple[dict[str, Any], dict[str, Path]]:
    config = read_json(config_path.resolve(strict=True))
    if (
        config.get("schema_version") != SCHEMA
        or config.get("source_dataset_id") != "gse296875"
        or config.get("target_dataset_id") != "gse244832"
        or config.get("models") != list(MODEL_IDS)
        or config.get("rotations") != list(ROTATIONS)
        or config.get("fit", {}).get("lsi_components") != 16
        or config.get("fit", {}).get("ridge_alpha") != 1.0
        or config.get("fit", {}).get("standardized_clip") != 8.0
    ):
        raise BlindBaselineRunError("campaign identity or fitted parameters differ")
    for field in (
        "condition_values_read",
        "phenotype_values_read",
        "rna_assay_read",
        "cross_assay_join_consumed",
        "sequence_features_consumed",
        "development_outcomes_read",
        "sealed_data_read",
        "metrics_calculated",
        "target_role_atac_consumed",
        "target_dataset_fit_or_adaptation",
        "missing_as_zero",
    ):
        if config.get("firewall", {}).get(field) is not False:
            raise BlindBaselineRunError(f"campaign firewall opened: {field}")
    if config.get("firewall", {}).get("prediction_commit_before_outcome_open") is not True:
        raise BlindBaselineRunError("prediction commit is not required")
    paths = {
        "registration": resolve_tree(
            root,
            config["execution_registration"],
            "gse244832_observed_atac_execution_registration",
        ),
        "axis": resolve_tree(
            root,
            config["exchange_axis"],
            "gse244832_reference_guarded_atac_exchange_axis",
        ),
    }
    for role in ("valid", "test"):
        paths[f"query_{role}"] = resolve_tree(
            root,
            config["query_authorities"][role],
            "gse244832_label_free_query_atac",
        )
        contract = read_json(paths[f"query_{role}"] / "query_contract.json")
        if (
            contract.get("context_role") != role
            or contract.get("scored_target_role") != ("test" if role == "valid" else "valid")
            or contract.get("fragment_cut_sites") != "start_plus_4_and_end_minus_5"
            or contract.get("query_transform_fit") is not False
            or contract.get("evaluator_outcome_artifact_read") is not False
        ):
            raise BlindBaselineRunError("query-ATAC contract differs")
    for name, record in config.get("implementation", {}).items():
        path = (root / str(record.get("path", ""))).resolve(strict=True)
        path.relative_to(root)
        if path.is_symlink() or not path.is_file() or digest(path) != record.get("sha256"):
            raise BlindBaselineRunError(f"implementation binding differs: {name}")
    return config, paths


def read_source_axis(path: Path) -> tuple[list[str], np.ndarray]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "donor_index",
            "donor_id",
            "lineage_index",
            "lineage_id",
        ):
            raise BlindBaselineRunError("source donor-lineage axis fields differ")
        rows = list(reader)
    if len(rows) != 156:
        raise BlindBaselineRunError("source donor-lineage axis size differs")
    donors: list[str] = []
    lineages: list[int] = []
    for row_index, row in enumerate(rows):
        donor_index, lineage_index = divmod(row_index, 4)
        if (
            int(row["donor_index"]) != donor_index
            or int(row["lineage_index"]) != lineage_index
            or row["lineage_id"] != LINEAGES[lineage_index]
        ):
            raise BlindBaselineRunError("source donor-lineage order differs")
        if lineage_index == 0:
            donors.append(row["donor_id"])
        lineages.append(lineage_index)
    if len(donors) != 39 or len(set(donors)) != 39:
        raise BlindBaselineRunError("source donor roster differs")
    return donors, np.asarray(lineages, dtype=np.int64)


def validate_source(
    source: Path,
    source_artifacts_sha256: str,
    config: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    manifest = verify_frozen_tree(source)
    contract = read_json(source / "contract.json")
    if (
        source.is_symlink()
        or digest(source / "ARTIFACTS.json") != source_artifacts_sha256
        or manifest.get("metadata", {}).get("artifact_class")
        != "gse296875_gse244832_source_window_counts"
        or contract.get("exchange_axis_artifacts_sha256")
        != config["exchange_axis"]["artifacts_sha256"]
        or contract.get("source_bigwigs_artifacts_sha256")
        != config["source_bigwigs"]["artifacts_sha256"]
        or contract.get("source_fragments_artifacts_sha256")
        != config["source_fragments"]["artifacts_sha256"]
        or contract.get("execution_registration_artifacts_sha256")
        != config["execution_registration"]["artifacts_sha256"]
        or contract.get("target_dataset_values_read") is not False
        or contract.get("metrics_calculated") is not False
        or contract.get("shape") != [39, 4, 2, 16000]
    ):
        raise BlindBaselineRunError("source fixed-window input differs")
    _donors, lineages = read_source_axis(source / "source_axis.tsv")
    counts = np.load(source / "source_counts.uint32.npy", mmap_mode="r", allow_pickle=False)
    if counts.shape != (39, 4, 2, 16_000) or counts.dtype != np.uint32:
        raise BlindBaselineRunError("source fixed-window count array differs")
    return counts, lineages


def read_query_axis(path: Path) -> np.ndarray:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "donor_index",
            "donor_id",
            "lineage_index",
            "lineage_id",
        ):
            raise BlindBaselineRunError("query donor-lineage axis fields differ")
        rows = list(reader)
    if len(rows) != 72:
        raise BlindBaselineRunError("query donor-lineage axis size differs")
    output = np.empty(72, dtype=np.int64)
    donors: list[str] = []
    for row_index, row in enumerate(rows):
        donor_index, lineage_index = divmod(row_index, 4)
        if (
            int(row["donor_index"]) != donor_index
            or int(row["lineage_index"]) != lineage_index
            or row["lineage_id"] != LINEAGES[lineage_index]
        ):
            raise BlindBaselineRunError("query donor-lineage order differs")
        output[row_index] = lineage_index
        if lineage_index == 0:
            donors.append(row["donor_id"])
    if len(donors) != 18 or len(set(donors)) != 18:
        raise BlindBaselineRunError("query donor roster differs")
    return output


def array_record(root: Path, path: Path, dtype: str) -> dict[str, str]:
    return {"path": str(path.relative_to(root)), "sha256": digest(path), "dtype": dtype}


def write_runtime_lock(output: Path) -> None:
    completed = subprocess.run(
        [sys.executable, "-m", "pip", "freeze"],
        check=True,
        capture_output=True,
        text=True,
    )
    write_text_exclusive(output / "pip_freeze.txt", completed.stdout)
    write_text_exclusive(
        output / "python_version.txt",
        f"{platform.python_version()}\n{sys.executable}\n",
    )


def write_prediction(
    *,
    output: Path,
    model_id: str,
    rotation_id: str,
    state: Any,
    predicted_eligible: np.ndarray,
    eligible: np.ndarray,
    unit_states: np.ndarray,
    config: Mapping[str, Any],
    paths: Mapping[str, Path],
    source: Path,
    source_artifacts_sha256: str,
) -> str:
    context_role, target_role = ROTATIONS[rotation_id]
    output.mkdir(parents=True, mode=0o750)
    state_path = output / "fit_state"
    export_state(state, state_path)
    predictions = np.full((18, 4, 16_000), np.nan, dtype=np.float32)
    flat_predictions = predictions.reshape(72, 16_000)
    flat_predictions[eligible] = np.asarray(predicted_eligible, dtype=np.float32)
    missing = np.broadcast_to(unit_states[:, :, None], predictions.shape).copy()
    predictions_path = output / "predictions.float32.npy"
    missing_path = output / "missing_state.uint8.npy"
    np.save(predictions_path, predictions, allow_pickle=False)
    np.save(missing_path, missing, allow_pickle=False)
    model_manifest = {
        "schema_version": "masld-bench-gse244832-task-native-model-execution-v1",
        "model_id": model_id,
        "rotation_id": rotation_id,
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse244832",
        "fit_scope": "gse296875_source_training_only",
        "source_training_target_counts_consumed": True,
        "source_input_artifacts_sha256": source_artifacts_sha256,
        "source_input_path": str(source),
        "execution_registration_artifacts_sha256": config["execution_registration"]["artifacts_sha256"],
        "exchange_axis_artifacts_sha256": config["exchange_axis"]["artifacts_sha256"],
        "context_role": context_role,
        "scored_target_role": target_role,
        "sequence_features_consumed": False,
        "rna_assay_consumed": False,
        "target_dataset_fit_or_adaptation": False,
        "per_donor_context_operation": "not_applicable",
        "query_transform_scope": "source_training_only",
        "mask_applied_before_query_transform": True,
        "ineligible_query_units_read_by_transform": False,
        "state_path": "fit_state",
        "state_npz_sha256": digest(state_path / "state.npz"),
        "state_metadata_sha256": digest(state_path / "metadata.json"),
        "metrics_calculated": False,
    }
    model_manifest_path = output / "model_execution_manifest.json"
    write_json_exclusive(model_manifest_path, model_manifest)
    bundle = {
        "schema_version": "masld-bench-atac-transport-prediction-bundle-v1",
        "bundle_id": f"{model_id}-source-fit-{rotation_id}-{source_artifacts_sha256[:16]}",
        "model_id": model_id,
        "input_regime": "observed_atac",
        "output_family": "masked_accessibility_count",
        "prediction_layout": "donor_lineage",
        "rotation_id": rotation_id,
        "axis_artifacts_sha256": config["exchange_axis"]["artifacts_sha256"],
        "execution_registration_artifacts_sha256": config["execution_registration"]["artifacts_sha256"],
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse244832",
        "query_atac_authority": {
            "path": str(paths[f"query_{context_role}"]),
            "artifacts_sha256": config["query_authorities"][context_role]["artifacts_sha256"],
            "context_role": context_role,
        },
        "model_manifest": {
            "path": "model_execution_manifest.json",
            "sha256": digest(model_manifest_path),
        },
        "arrays": {
            "predictions": array_record(output, predictions_path, "float32"),
            "missing_state": array_record(output, missing_path, "uint8"),
        },
        "prediction_scale": "source_fitted_expected_deduplicated_tn5_insertion_count",
        "seeds": [20260825],
        "evidence_contract": {
            "condition_labels_consumed": False,
            "phenotype_values_consumed": False,
            "rna_assay_consumed": False,
            "cross_assay_join_consumed": False,
            "sequence_features_consumed": False,
            "source_native_reference_bundle_assumed_resolved": False,
            "development_outcomes_consumed": False,
            "sealed_data_consumed": False,
            "target_role_atac_consumed": False,
            "target_role_contigs_entered_query_transform": False,
            "missing_as_zero": False,
            "scored_target_role": target_role,
            "observed_atac_context_role": context_role,
            "query_atac_consumed": True,
            "receptive_field_bp": "not_applicable_global_observed_atac_encoder",
            "query_transform_scope": "source_training_only",
            "per_donor_context_operation": "not_applicable",
        },
        "metrics_calculated": False,
        "champion_claim_allowed": False,
    }
    write_json_exclusive(output / "prediction_bundle.json", bundle)
    write_json_exclusive(
        output / "fit_predict_receipt.json",
        {
            "schema_version": "masld-bench-gse244832-task-native-fit-predict-receipt-v1",
            "status": "passed",
            "model_id": model_id,
            "rotation_id": rotation_id,
            "source_donor_lineage_units": 156,
            "query_donor_lineage_units": 72,
            "eligible_query_units_transformed": int(eligible.sum()),
            "ineligible_query_units_transformed": 0,
            "mask_applied_before_query_transform": True,
            "target_role_atac_consumed": False,
            "condition_or_phenotype_read": False,
            "development_outcomes_read": False,
            "metrics_calculated": False,
        },
    )
    return freeze_tree(
        output,
        {
            "artifact_class": "gse244832_task_native_blind_prediction",
            "model_id": model_id,
            "rotation_id": rotation_id,
            "source_dataset_id": "gse296875",
            "target_dataset_id": "gse244832",
            "target_role_atac_consumed": False,
            "metrics_calculated": False,
            "status": "passed",
        },
    )


def run(
    root: Path,
    config_path: Path,
    source_input: Path,
    source_artifacts_sha256: str,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise BlindBaselineRunError("refusing to overwrite blind baseline run")
    root = root.resolve(strict=True)
    config, paths = validate_config(root, config_path)
    source = source_input.resolve(strict=True)
    source.relative_to(root / "executions")
    source_counts, source_lineages = validate_source(
        source, source_artifacts_sha256, config
    )
    registration = read_json(paths["registration"] / "execution_contract.json")
    unit_states = np.load(
        paths["registration"] / registration["query_lineage_missing_state_path"],
        allow_pickle=False,
    )
    if (
        unit_states.shape != (18, 4)
        or unit_states.dtype != np.uint8
        or set(np.unique(unit_states)) != {0, 1, 3}
        or int((unit_states == 0).sum()) != 48
        or int((unit_states == 1).sum()) != 4
        or int((unit_states == 3).sum()) != 20
    ):
        raise BlindBaselineRunError("registered query missingness mask differs")
    eligible = unit_states.reshape(-1) == 0
    output.mkdir(parents=True, mode=0o750)
    prediction_records: list[dict[str, Any]] = []
    commit_records: list[dict[str, Any]] = []
    for rotation_id, (context_role, target_role) in ROTATIONS.items():
        query_root = paths[f"query_{context_role}"]
        query_lineages = read_query_axis(query_root / "query_axis.tsv")
        query_memmap = np.load(
            query_root / "query_unit_fragment_total.uint32.npy",
            mmap_mode="r",
            allow_pickle=False,
        )
        if query_memmap.shape != (18, 4, 16_000) or query_memmap.dtype != np.uint32:
            raise BlindBaselineRunError("query fixed-window count array differs")
        query_context = np.asarray(
            query_memmap.reshape(72, 16_000)[eligible], dtype=np.float64
        )
        eligible_lineages = query_lineages[eligible]
        source_context = np.asarray(
            source_counts[:, :, ROLE_INDEX[context_role], :], dtype=np.float64
        ).reshape(156, 16_000)
        source_target = np.asarray(
            source_counts[:, :, ROLE_INDEX[target_role], :], dtype=np.float64
        ).reshape(156, 16_000)
        for model_id in MODEL_IDS:
            state = fit(
                model_id,
                context_counts=source_context,
                target_counts=source_target,
                lineage_indices=source_lineages,
                lsi_components=int(config["fit"]["lsi_components"]),
                ridge_alpha=float(config["fit"]["ridge_alpha"]),
                standardized_clip=float(config["fit"]["standardized_clip"]),
            )
            predicted = predict(
                state,
                context_counts=query_context,
                lineage_indices=eligible_lineages,
            )
            prediction_root = output / "predictions" / model_id / rotation_id
            prediction_sha = write_prediction(
                output=prediction_root,
                model_id=model_id,
                rotation_id=rotation_id,
                state=state,
                predicted_eligible=predicted,
                eligible=eligible,
                unit_states=unit_states,
                config=config,
                paths=paths,
                source=source,
                source_artifacts_sha256=source_artifacts_sha256,
            )
            prediction_records.append(
                {
                    "model_id": model_id,
                    "rotation_id": rotation_id,
                    "path": str(prediction_root.relative_to(output)),
                    "artifacts_sha256": prediction_sha,
                }
            )
            commit_root = output / "commits" / model_id / rotation_id
            commit_receipt = commit(
                prediction_root=prediction_root,
                prediction_artifacts_sha256=prediction_sha,
                axis_root=paths["axis"],
                axis_artifacts_sha256=config["exchange_axis"]["artifacts_sha256"],
                registration_root=paths["registration"],
                registration_artifacts_sha256=config["execution_registration"]["artifacts_sha256"],
                output=commit_root,
            )
            commit_records.append(
                {
                    "model_id": model_id,
                    "rotation_id": rotation_id,
                    "path": str(commit_root.relative_to(output)),
                    "artifacts_sha256": digest(commit_root / "ARTIFACTS.json"),
                    "status": commit_receipt["status"],
                }
            )
    write_json_exclusive(
        output / "receipt.json",
        {
            "schema_version": "masld-bench-gse244832-task-native-blind-baseline-run-receipt-v1",
            "status": "passed",
            "source_artifacts_sha256": source_artifacts_sha256,
            "execution_registration_artifacts_sha256": config["execution_registration"]["artifacts_sha256"],
            "prediction_records": prediction_records,
            "commit_records": commit_records,
            "predictions": len(prediction_records),
            "commits": len(commit_records),
            "all_predictions_committed_before_outcome_open": True,
            "eligible_query_units": 48,
            "structurally_missing_query_units": 4,
            "below_qc_query_units": 20,
            "condition_or_phenotype_read": False,
            "target_role_atac_consumed": False,
            "development_outcomes_read": False,
            "sealed_data_read": False,
            "metrics_calculated": False,
            "champion_claim_allowed": False,
        },
    )
    write_runtime_lock(output)
    artifact_sha = freeze_tree(
        output,
        {
            "artifact_class": "gse244832_task_native_blind_baseline_run",
            "models": list(MODEL_IDS),
            "rotations": list(ROTATIONS),
            "prediction_commits": 6,
            "development_outcomes_read": False,
            "metrics_calculated": False,
            "status": "passed",
        },
    )
    return {"output": str(output), "artifacts_sha256": artifact_sha, "commits": 6}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", dest="config_path", required=True, type=Path)
    parser.add_argument("--source-input", required=True, type=Path)
    parser.add_argument("--source-artifacts-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    print(json.dumps(run(**vars(arguments)), sort_keys=True))


if __name__ == "__main__":
    main()
